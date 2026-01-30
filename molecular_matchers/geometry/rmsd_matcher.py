#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD-based 3D geometry molecular matcher.

Compares molecular structures using Root Mean Square Deviation (RMSD)
after optimal alignment (Kabsch algorithm).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..base import MoleculeMatcher
from ..registry import register_matcher


@register_matcher("rmsd")
class RMSDMatcher(MoleculeMatcher):
    """
    RMSD-based matcher for 3D structural comparison.

    Config:
        threshold: float = 0.5        # RMSD threshold in Angstrom
        allow_reorder: bool = True    # Allow atom reordering (Hungarian algorithm)
        heavy_only: bool = False      # Compare only heavy atoms (ignore H)

    Uses Kabsch algorithm for optimal alignment and computes RMSD.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self._threshold = self.config.get("threshold", 0.5)  # Angstrom
        self._allow_reorder = self.config.get("allow_reorder", False)
        self._heavy_only = self.config.get("heavy_only", False)

    @property
    def name(self) -> str:
        return "rmsd"

    def _read_xyz(self, xyz_file: Path) -> Optional[Tuple[List[str], np.ndarray]]:
        """
        Read XYZ file and extract elements and coordinates.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Tuple of (elements list, coordinates array Nx3), or None if read fails
        """
        try:
            with open(xyz_file) as f:
                lines = f.readlines()

            natoms = int(lines[0].strip())
            elements = []
            coords = []

            for i in range(2, 2 + natoms):
                parts = lines[i].split()
                elements.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

            return elements, np.array(coords)

        except Exception:
            return None

    def _filter_heavy_atoms(
        self, elements: List[str], coords: np.ndarray
    ) -> Tuple[List[str], np.ndarray]:
        """Filter to keep only heavy atoms (non-hydrogen)."""
        heavy_indices = [i for i, el in enumerate(elements) if el != "H"]
        heavy_elements = [elements[i] for i in heavy_indices]
        heavy_coords = coords[heavy_indices]
        return heavy_elements, heavy_coords

    def _kabsch_align(self, coords_ref: np.ndarray, coords_mob: np.ndarray) -> np.ndarray:
        """
        Align coords_mob to coords_ref using Kabsch algorithm.

        Args:
            coords_ref: Reference coordinates (Nx3)
            coords_mob: Mobile coordinates to be aligned (Nx3)

        Returns:
            Aligned mobile coordinates (Nx3)
        """
        # Center both coordinate sets
        centroid_ref = np.mean(coords_ref, axis=0)
        centroid_mob = np.mean(coords_mob, axis=0)

        coords_ref_centered = coords_ref - centroid_ref
        coords_mob_centered = coords_mob - centroid_mob

        # Compute covariance matrix
        H = coords_mob_centered.T @ coords_ref_centered

        # SVD
        U, S, Vt = np.linalg.svd(H)

        # Compute rotation matrix
        R = Vt.T @ U.T

        # Handle reflection case
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        # Rotate and translate
        aligned = coords_mob_centered @ R.T + centroid_ref
        return aligned

    def _compute_rmsd(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """
        Compute RMSD between two aligned coordinate sets.

        Args:
            coords1: First coordinate set (Nx3)
            coords2: Second coordinate set (Nx3)

        Returns:
            RMSD value in Angstrom
        """
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    def _hungarian_match(
        self, elements1: List[str], coords1: np.ndarray, elements2: List[str], coords2: np.ndarray
    ) -> np.ndarray:
        """
        Find optimal atom correspondence using Hungarian algorithm.

        Args:
            elements1: Elements of first structure
            coords1: Coordinates of first structure
            elements2: Elements of second structure
            coords2: Coordinates of second structure

        Returns:
            Permutation indices to reorder coords2 to match coords1
        """
        try:
            from scipy.optimize import linear_sum_assignment
        except ImportError:
            # Fallback: simple index-based matching
            return np.arange(len(elements1))

        # Build cost matrix (distance between atoms of same element)
        n1 = len(elements1)
        n2 = len(elements2)
        cost_matrix = np.full((n1, n2), np.inf)

        for i, el1 in enumerate(elements1):
            for j, el2 in enumerate(elements2):
                if el1 == el2:
                    cost_matrix[i, j] = np.linalg.norm(coords1[i] - coords2[j])

        # Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        return col_ind

    def compute_signature(self, xyz_file: Path) -> Optional[Any]:
        """
        Read XYZ file and return coordinates as signature.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Tuple of (elements, coordinates), or None if read fails
        """
        result = self._read_xyz(xyz_file)
        if result is None:
            return None

        elements, coords = result

        if self._heavy_only:
            elements, coords = self._filter_heavy_atoms(elements, coords)

        return (elements, coords)

    def compare(self, sig1: Any, sig2: Any) -> bool:
        """
        Compare two structures using RMSD after optimal alignment.

        Args:
            sig1: First signature (elements, coords)
            sig2: Second signature (elements, coords)

        Returns:
            True if RMSD < threshold
        """
        if sig1 is None or sig2 is None:
            return False

        try:
            elements1, coords1 = sig1
            elements2, coords2 = sig2

            # Check atom count
            if len(elements1) != len(elements2):
                return False

            # Optionally reorder atoms for best match
            if self._allow_reorder:
                perm = self._hungarian_match(elements1, coords1, elements2, coords2)
                coords2 = coords2[perm]
                elements2 = [elements2[i] for i in perm]

            # Check element composition
            if set(elements1) != set(elements2):
                return False

            # Align and compute RMSD
            coords2_aligned = self._kabsch_align(coords1, coords2)
            rmsd = self._compute_rmsd(coords1, coords2_aligned)

            return rmsd < self._threshold

        except Exception:
            return False

    def _serialize_signature(self, sig: Any) -> Any:
        """Serialize signature for JSON output."""
        if sig is None:
            return None
        elements, coords = sig
        return {
            "elements": elements,
            "natoms": len(elements),
            "coords_shape": coords.shape,
        }

    def cross_isomorphism(
        self,
        true_r: Path,
        true_p: Path,
        irc_r: Path,
        irc_p: Path,
    ) -> Any:
        """
        Cross-compare using RMSD with computed RMSD values in metadata.

        Args:
            true_r: Original reactant XYZ file
            true_p: Original product XYZ file
            irc_r: IRC backward endpoint XYZ file
            irc_p: IRC forward endpoint XYZ file

        Returns:
            MatchResult with RMSD values in metadata
        """
        from ..base import MatchResult

        try:
            sig_true_r = self.compute_signature(true_r)
            sig_true_p = self.compute_signature(true_p)
            sig_irc_r = self.compute_signature(irc_r)
            sig_irc_p = self.compute_signature(irc_p)

            if None in (sig_true_r, sig_true_p, sig_irc_r, sig_irc_p):
                return MatchResult(
                    performed=True,
                    success=False,
                    endpoint_match="Error",
                    method=self.name,
                    error="Failed to read one or more XYZ files",
                )

            # Compute RMSD values for metadata
            def compute_rmsd_value(sig1, sig2):
                """Helper to compute actual RMSD value."""
                try:
                    _, coords1 = sig1
                    _, coords2 = sig2
                    if self._allow_reorder:
                        elements1, coords1 = sig1
                        elements2, coords2 = sig2
                        perm = self._hungarian_match(elements1, coords1, elements2, coords2)
                        coords2 = coords2[perm]
                    coords2_aligned = self._kabsch_align(coords1, coords2)
                    return self._compute_rmsd(coords1, coords2_aligned)
                except Exception:
                    return None

            rmsd_true_r_irc_r = compute_rmsd_value(sig_true_r, sig_irc_r)
            rmsd_true_p_irc_p = compute_rmsd_value(sig_true_p, sig_irc_p)
            rmsd_true_r_irc_p = compute_rmsd_value(sig_true_r, sig_irc_p)
            rmsd_true_p_irc_r = compute_rmsd_value(sig_true_p, sig_irc_r)
            rmsd_irc_r_irc_p = compute_rmsd_value(sig_irc_r, sig_irc_p)

            # Cross-compare
            matches = (
                self.compare(sig_true_r, sig_irc_r),
                self.compare(sig_true_p, sig_irc_p),
                self.compare(sig_true_r, sig_irc_p),
                self.compare(sig_true_p, sig_irc_r),
            )

            # Reaction type
            rxn_status = (
                "Conformational change"
                if self.compare(sig_irc_r, sig_irc_p)
                else "Chemical reaction"
            )

            # Endpoint match category
            if matches in [(True, True, False, False), (False, False, True, True)]:
                endpoint_match = "2-end match"
            elif any(matches):
                endpoint_match = "1-end match"
            else:
                endpoint_match = "No match"

            return MatchResult(
                performed=True,
                success=True,
                endpoint_match=endpoint_match,
                matches=tuple(bool(m) for m in matches),
                rxn_status=rxn_status,
                method=self.name,
                metadata={
                    "threshold_angstrom": float(self._threshold),
                    "rmsd_true_r_irc_r": float(rmsd_true_r_irc_r) if rmsd_true_r_irc_r is not None else None,
                    "rmsd_true_p_irc_p": float(rmsd_true_p_irc_p) if rmsd_true_p_irc_p is not None else None,
                    "rmsd_true_r_irc_p": float(rmsd_true_r_irc_p) if rmsd_true_r_irc_p is not None else None,
                    "rmsd_true_p_irc_r": float(rmsd_true_p_irc_r) if rmsd_true_p_irc_r is not None else None,
                    "rmsd_irc_r_irc_p": float(rmsd_irc_r_irc_p) if rmsd_irc_r_irc_p is not None else None,
                    "heavy_only": bool(self._heavy_only),
                    "allow_reorder": bool(self._allow_reorder),
                },
            )

        except Exception as e:
            return MatchResult(
                performed=True,
                success=False,
                endpoint_match="Error",
                method=self.name,
                error=str(e),
            )
