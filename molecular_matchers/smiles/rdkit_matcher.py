#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RDKit SMILES-based molecular matcher.

Converts XYZ to canonical SMILES using RDKit with multiple fallback strategies
(Hueckel, Vdw, Basic) for bond determination.

Supports optional atom map numbers for tracking atoms across reactions.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from ..base import MoleculeMatcher
from ..registry import register_matcher


@register_matcher("smiles_rdkit")
class RDKitSmilesMatcher(MoleculeMatcher):
    """
    SMILES-based matcher using RDKit with robust bond determination.

    Config:
        atom_map: bool = False  # Include atom map numbers in SMILES (AAM)

    Uses multiple fallback strategies for bond determination:
    - Hueckel: Use Hueckel orbital theory
    - Vdw: Use van der Waals radii
    - Basic: Basic distance-based

    Always generates canonical SMILES with stereochemistry.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self._atom_map = self.config.get("atom_map", False)
        self._all_hs_explicit = self.config.get("all_hs_explicit", True)

    @property
    def name(self) -> str:
        return "smiles_rdkit"

    def _check_dependencies(self) -> bool:
        """Check if RDKit is available."""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
            from rdkit.Chem import rdmolops
            return True
        except Exception:
            return False

    def _check_charge(self, smiles_text: str) -> bool:
        """
        Check if SMILES represents a neutral molecule.

        Args:
            smiles_text: SMILES string (may contain dots for multi-component)

        Returns:
            True if all components are neutral
        """
        from rdkit import Chem

        for part in smiles_text.split('.'):
            mol = Chem.MolFromSmiles(part)
            if mol is None:
                return False
            if sum(atom.GetFormalCharge() for atom in mol.GetAtoms()) != 0:
                return False
        return True

    def _post_process(self, smiles_text: str, xyz_path: Path, use_hueckel: bool) -> str:
        """
        Post-process SMILES to fix charge issues using original XYZ coordinates.

        Args:
            smiles_text: Initial SMILES from RDKit
            xyz_path: Path to original XYZ file
            use_hueckel: Whether Hueckel strategy was used

        Returns:
            Corrected SMILES string

        Raises:
            ValueError: If post-processing fails
        """
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds

        # Simple XYZ parser (avoid ASE dependency)
        class SimpleStructure:
            def __init__(self, atoms, positions):
                self.atoms = atoms
                self._positions = positions

            def get_positions(self):
                return self._positions

        def read(xyz_path):
            with open(xyz_path) as f:
                lines = f.readlines()
            natoms = int(lines[0].strip())
            atoms = []
            positions = []
            for i in range(2, 2 + natoms):
                parts = lines[i].split()
                atoms.append(parts[0])
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
            return SimpleStructure(atoms, positions)

        params = Chem.SmilesParserParams()
        params.removeHs = False

        if self._check_charge(smiles_text):
            return smiles_text

        structure = read(xyz_path)
        pos = structure.get_positions()
        out_answer = []

        for part in smiles_text.split('.'):
            mol = Chem.MolFromSmiles(part, params)
            if mol is None:
                raise ValueError("SMILES parse failed")

            am_list = []
            xyz_lines = [f"{mol.GetNumAtoms()}\n"]
            for atom in mol.GetAtoms():
                this_am = atom.GetAtomMapNum()
                am_list.append(this_am)
                if this_am <= 0 or this_am > len(pos):
                    raise ValueError("Atom map number out of range")
                x, y, z = pos[this_am - 1]
                xyz_lines.append(f"{atom.GetSymbol()} {x} {y} {z}")

            xyz_block = '\n'.join(xyz_lines)

            # Build strategies for post-processing
            strategies = []
            if use_hueckel:
                strategies.append({"params": {"useHueckel": True, "useAtomMap": False}, "name": "Hueckel"})
            strategies.extend([
                {"params": {"useVdw": True, "useAtomMap": False}, "name": "Vdw"},
                {"params": {"useAtomMap": False}, "name": "Basic"}
            ])
            if not use_hueckel:
                strategies.append({"params": {"useHueckel": True, "useAtomMap": False}, "name": "Hueckel"})

            last_error = None
            mini_mol = None

            for strategy in strategies:
                try:
                    attempt = Chem.MolFromXYZBlock(xyz_block)
                    if attempt is None:
                        raise ValueError("XYZ block parse failed")
                    rdDetermineBonds.DetermineBonds(attempt, **strategy["params"])
                    mini_mol = attempt
                    break
                except Exception as e:
                    last_error = e
                    continue

            if mini_mol is None:
                raise ValueError(f"DetermineBonds failed: {last_error}")

            for atom in mini_mol.GetAtoms():
                atom.SetAtomMapNum(am_list[atom.GetIdx()])
            out_answer.append(Chem.MolToSmiles(mini_mol, canonical=True, isomericSmiles=True))

        return '.'.join(out_answer)

    def compute_signature(self, xyz_file: Path) -> Optional[str]:
        """
        Convert XYZ to canonical SMILES using RDKit with multiple fallback strategies.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Canonical SMILES string, or None if conversion fails
        """
        if not self._check_dependencies():
            return None

        if not os.path.exists(xyz_file):
            return None

        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
            from rdkit.Chem import rdmolops

            # Step 1: Read raw molecule from XYZ
            raw_mol = Chem.MolFromXYZFile(str(xyz_file))
            if raw_mol is None:
                return None

            # Define multiple strategies
            strategies = [
                {"params": {"useHueckel": True, "useAtomMap": True}, "name": "Hueckel"},
                {"params": {"useVdw": True, "useAtomMap": True}, "name": "Vdw"},
                {"params": {"useAtomMap": True}, "name": "Basic"}
            ]

            mol = None
            saved_map_nums = []
            success = False
            last_error = None
            use_hueckel = False

            # Step 2: Try each strategy in order
            for strategy in strategies:
                try:
                    mol = Chem.Mol(raw_mol)
                    rdDetermineBonds.DetermineBonds(mol, **strategy["params"])

                    # Backup and clear MapNums (for proper Sanitization)
                    current_saved_maps = []
                    for atom in mol.GetAtoms():
                        current_saved_maps.append(atom.GetAtomMapNum())
                        atom.SetAtomMapNum(0)

                    # Step 3: Sanitize
                    Chem.SanitizeMol(mol)

                    saved_map_nums = current_saved_maps
                    success = True
                    use_hueckel = bool(strategy["params"].get("useHueckel", False))
                    break

                except Exception as e:
                    last_error = e
                    continue

            if not success:
                return None

            # Step 4: Assign stereochemistry from 3D coordinates
            rdmolops.AssignStereochemistryFrom3D(mol)

            # Step 5: Restore AtomMapNum if atom_map is enabled
            if self._atom_map:
                for i, atom in enumerate(mol.GetAtoms()):
                    atom.SetAtomMapNum(saved_map_nums[i])

            # Step 6: Generate SMILES (always canonical with stereo)
            # Remove Hs if all_hs_explicit=False
            mol_for_smiles = mol
            if not self._all_hs_explicit:
                mol_for_smiles = Chem.RemoveHs(mol)
                if self._atom_map:
                    # Re-assign map nums after RemoveHs (indices change)
                    for atom in mol_for_smiles.GetAtoms():
                        atom.SetAtomMapNum(atom.GetIdx() + 1)

            aam_smi = Chem.MolToSmiles(mol_for_smiles, canonical=True, isomericSmiles=True)

            # Step 7: Post-process to fix charge issues
            try:
                processed_aam = self._post_process(aam_smi, xyz_file, use_hueckel)
            except Exception:
                # Fall back to unprocessed SMILES if post-process fails
                processed_aam = aam_smi

            return processed_aam

        except Exception:
            return None

    def compare(self, sig1: str, sig2: str) -> bool:
        """
        Compare two SMILES strings for equality.

        Args:
            sig1: First SMILES string
            sig2: Second SMILES string

        Returns:
            True if SMILES are identical
        """
        if sig1 is None or sig2 is None:
            return False
        return sig1 == sig2

    def _build_cross_isomorphism_metadata(
        self, sig_true_r: Any, sig_true_p: Any, sig_irc_r: Any, sig_irc_p: Any
    ) -> Dict[str, Any]:
        """Add atom_map and all_hs_explicit settings to metadata along with signatures."""
        metadata = super()._build_cross_isomorphism_metadata(sig_true_r, sig_true_p, sig_irc_r, sig_irc_p)
        metadata["atom_map"] = self._atom_map
        metadata["all_hs_explicit"] = self._all_hs_explicit
        return metadata
