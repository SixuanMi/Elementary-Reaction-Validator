#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenBabel SMILES-based molecular matcher.

Converts XYZ to canonical SMILES (including stereochemistry) using OpenBabel/pybel.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from ..base import MoleculeMatcher
from ..registry import register_matcher


@register_matcher("smiles_openbabel")
class OpenBabelSmilesMatcher(MoleculeMatcher):
    """
    SMILES-based matcher using OpenBabel/pybel.

    Returns canonical SMILES string (with stereochemistry) as signature.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)

    @property
    def name(self) -> str:
        return "smiles_openbabel"

    def _check_dependencies(self) -> bool:
        """Check if OpenBabel/pybel is available."""
        try:
            from openbabel import pybel  # type: ignore
            return True
        except Exception:
            return False

    def compute_signature(self, xyz_file: Path) -> Optional[str]:
        """
        Convert XYZ to canonical SMILES using OpenBabel/pybel.

        Args:
            xyz_file: Path to XYZ file

        Returns:
            Canonical SMILES string with stereochemistry, or None if conversion fails
        """
        if not self._check_dependencies():
            return None

        try:
            from openbabel import pybel  # type: ignore
            mol = next(pybel.readfile("xyz", str(xyz_file)))
            return mol.write("can").strip().split()[0]
        except Exception:
            return None

    def compare(self, sig1: str, sig2: str) -> bool:
        """
        Compare two SMILES strings for equality.

        Args:
            sig1: First SMILES string
            sig2: Second SMILES string

        Returns:
            True if SMILES are identical (canonical comparison)
        """
        if sig1 is None or sig2 is None:
            return False
        return sig1 == sig2
