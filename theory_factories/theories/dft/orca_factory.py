#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ORCA theory factory.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from ...base import TheoryFactoryBase
from ...registry import register_theory


@register_theory("orca")
class ORCAFactory(TheoryFactoryBase):
    """
    ORCA theory factory.

    Config:
        {
            "orcasimpleinput": "BP86 def2-SVP",  # ORCA simple input line
            "orcablocks": "",                    # ORCA blocks (optional)
            "moreadfile": None                   # Path to .molden file (optional)
        }
    """

    # Hardcoded ORCA directory for parallel runs
    _ORCADIR = "/inspire/hdd/project/chemicalreaction/misixuan-CZXS24220243/soft/orca_6_0_1_linux_x86-64_shared_openmpi416_avx2"

    @property
    def name(self) -> str:
        return "orca"

    def validate_config(self) -> bool:
        """Validate ORCA configuration."""
        return "orcasimpleinput" in self.config

    def create_theory(
        self, numcores: int = 1, printlevel: int = 1, **kwargs
    ) -> Any:
        """
        Create ORCA theory instance.

        Args:
            numcores: Number of CPU cores
            printlevel: Print level for output

        Returns:
            ORCATheory instance
        """
        from ash import ORCATheory

        if not self.validate_config():
            raise ValueError(f"Invalid ORCA config: {self.config}")

        return ORCATheory(
            orcadir=self._ORCADIR,
            orcasimpleinput=self.config["orcasimpleinput"],
            orcablocks=self.config.get("orcablocks", ""),
            numcores=numcores,
            printlevel=printlevel,
            filename="orca_",
            moreadfile=self.config.get("moreadfile"),
        )

    def get_default_config(self) -> Dict[str, Any]:
        """Return default ORCA configuration."""
        return {
            "orcasimpleinput": "BP86 def2-SVP",
            "orcablocks": "",
        }
