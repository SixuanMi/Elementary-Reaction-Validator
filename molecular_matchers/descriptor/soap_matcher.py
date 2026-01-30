#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SOAP (Smooth Overlap of Atomic Positions) descriptor-based molecular matcher.

Compares molecular structures using SOAP descriptors from dscribe with
AverageKernel for similarity comparison.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np

from ..base import MoleculeMatcher
from ..registry import register_matcher


@register_matcher("soap")
class SOAPMatcher(MoleculeMatcher):
    """
    SOAP descriptor-based matcher using dscribe AverageKernel.

    Config:
        r_cut: float = 6.0                # Cutoff radius in Angstrom
        n_max: int = 8                    # Number of radial basis functions
        l_max: int = 6                    # Maximum angular momentum
        sigma: float = 0.1                # Width of Gaussian (smoothness)
        periodic: bool = False            # Periodic boundary conditions
        kernel_metric: str = "linear"     # "linear" or "rbf" for AverageKernel
        kernel_gamma: float = None        # RBF gamma (None=auto, only for rbf metric)
        threshold_similarity: float = 0.99  # Similarity threshold for match

    Uses dscribe SOAP descriptor with AverageKernel for structure comparison.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self._r_cut = self.config.get("r_cut", 6.0)
        self._n_max = self.config.get("n_max", 8)
        self._l_max = self.config.get("l_max", 6)
        self._sigma = self.config.get("sigma", 0.1)
        self._periodic = self.config.get("periodic", False)
        self._kernel_metric = self.config.get("kernel_metric", "linear")
        self._kernel_gamma = self.config.get("kernel_gamma", None)
        self._threshold = self.config.get("threshold_similarity", 0.99)

    @property
    def name(self) -> str:
        return "soap"

    def _check_dependencies(self) -> bool:
        """Check if dscribe is available."""
        try:
            from dscribe.descriptors import SOAP
            from dscribe.kernels import AverageKernel
            return True
        except Exception:
            return False

    def _read_xyz_to_ase(self, xyz_file: Path):
        """
        Read XYZ file and create ASE Atoms object.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            ASE Atoms object, or None if read fails
        """
        try:
            from ase import Atoms

            with open(xyz_file) as f:
                lines = f.readlines()

            natoms = int(lines[0].strip())
            symbols = []
            positions = []

            for i in range(2, 2 + natoms):
                parts = lines[i].split()
                symbols.append(parts[0])
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])

            return Atoms(symbols=symbols, positions=positions)

        except Exception:
            return None

    def _get_soap_descriptor_and_kernel(self, xyz_files: list) -> Tuple[Any, Any, list]:
        """
        Create SOAP descriptor and AverageKernel for the given structures.

        Args:
            xyz_files: List of XYZ file paths

        Returns:
            Tuple of (SOAP object, AverageKernel object, list of descriptors)
        """
        if not self._check_dependencies():
            return None, None, []

        try:
            from dscribe.descriptors import SOAP
            from dscribe.kernels import AverageKernel
            from ase import Atoms

            # Read all structures to get unique species
            all_species = set()
            structures = []
            for xyz_file in xyz_files:
                atoms = self._read_xyz_to_ase(xyz_file)
                if atoms is None:
                    return None, None, []
                structures.append(atoms)
                all_species.update(atoms.get_chemical_symbols())

            species_list = sorted(list(all_species))

            # Create SOAP descriptor
            soap = SOAP(
                species=species_list,
                r_cut=self._r_cut,
                n_max=self._n_max,
                l_max=self._l_max,
                sigma=self._sigma,
                periodic=self._periodic,
                compression={"mode": "off"},
                sparse=False,
            )

            # Create descriptors for all structures
            descriptors = []
            for struct in structures:
                desc = soap.create(struct)
                descriptors.append(desc)

            # Create AverageKernel
            kernel_params = {"metric": self._kernel_metric}
            if self._kernel_metric == "rbf" and self._kernel_gamma is not None:
                kernel_params["gamma"] = self._kernel_gamma

            kernel = AverageKernel(**kernel_params)

            return soap, kernel, descriptors

        except Exception as e:
            return None, None, []

    def compute_signature(self, xyz_file: Path) -> Optional[Any]:
        """
        Read XYZ file for SOAP descriptor computation.

        Note: Actual SOAP computation is done in cross_isomorphism for efficiency.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Path to XYZ file (for later batch processing)
        """
        if not self._check_dependencies():
            return None

        if not xyz_file.exists():
            return None

        # Just store the path for batch processing
        return xyz_file

    def compare(self, sig1: Any, sig2: Any) -> bool:
        """
        Compare two structures - not used for SOAP (batch processing in cross_isomorphism).

        Args:
            sig1: First signature (XYZ path)
            sig2: Second signature (XYZ path)

        Returns:
            True if similarity >= threshold
        """
        # For individual comparison, we need to compute descriptors
        if sig1 is None or sig2 is None:
            return False

        try:
            soap, kernel, descriptors = self._get_soap_descriptor_and_kernel([sig1, sig2])
            if descriptors is None or len(descriptors) < 2:
                return False

            kernel_matrix = kernel.create(descriptors)
            similarity = float(kernel_matrix[0, 1])

            return similarity >= self._threshold

        except Exception:
            return False

    def _serialize_signature(self, sig: Any) -> Any:
        """Serialize signature for JSON output."""
        if sig is None:
            return None
        return {"xyz_path": str(sig)}

    def cross_isomorphism(
        self,
        true_r: Path,
        true_p: Path,
        irc_r: Path,
        irc_p: Path,
    ) -> Any:
        """
        Cross-compare using SOAP with similarity values in metadata.

        Uses dscribe AverageKernel for efficient batch computation.

        Args:
            true_r: Original reactant XYZ file
            true_p: Original product XYZ file
            irc_r: IRC backward endpoint XYZ file
            irc_p: IRC forward endpoint XYZ file

        Returns:
            MatchResult with similarity values in metadata
        """
        from ..base import MatchResult

        try:
            # Batch compute SOAP descriptors and kernel matrix
            xyz_files = [true_r, true_p, irc_r, irc_p]
            soap, kernel, descriptors = self._get_soap_descriptor_and_kernel(xyz_files)

            if descriptors is None or len(descriptors) < 4:
                return MatchResult(
                    performed=True,
                    success=False,
                    endpoint_match="Error",
                    method=self.name,
                    error="Failed to compute SOAP descriptors",
                )

            # Compute kernel matrix
            kernel_matrix = kernel.create(descriptors)

            # Extract similarities (matrix is symmetric)
            # Order: true_r, true_p, irc_r, irc_p
            # Indices:    0       1       2       3
            sim_true_r_irc_r = float(kernel_matrix[0, 2])
            sim_true_p_irc_p = float(kernel_matrix[1, 3])
            sim_true_r_irc_p = float(kernel_matrix[0, 3])
            sim_true_p_irc_r = float(kernel_matrix[1, 2])
            sim_irc_r_irc_p = float(kernel_matrix[2, 3])

            # Cross-compare
            matches = (
                sim_true_r_irc_r >= self._threshold,
                sim_true_p_irc_p >= self._threshold,
                sim_true_r_irc_p >= self._threshold,
                sim_true_p_irc_r >= self._threshold,
            )

            # Reaction type
            rxn_status = (
                "Conformational change"
                if sim_irc_r_irc_p >= self._threshold
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
                matches=matches,
                rxn_status=rxn_status,
                method=self.name,
                metadata={
                    "threshold": self._threshold,
                    "r_cut": self._r_cut,
                    "n_max": self._n_max,
                    "l_max": self._l_max,
                    "sigma": self._sigma,
                    "kernel_metric": self._kernel_metric,
                    "kernel_gamma": self._kernel_gamma,
                    "similarity_true_r_irc_r": sim_true_r_irc_r,
                    "similarity_true_p_irc_p": sim_true_p_irc_p,
                    "similarity_true_r_irc_p": sim_true_r_irc_p,
                    "similarity_true_p_irc_r": sim_true_p_irc_r,
                    "similarity_irc_r_irc_p": sim_irc_r_irc_p,
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
