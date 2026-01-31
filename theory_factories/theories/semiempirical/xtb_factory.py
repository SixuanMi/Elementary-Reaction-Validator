#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xTB theory factory.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from ...base import TheoryFactoryBase
from ...registry import register_theory


@register_theory("xtb")
class XTBFactory(TheoryFactoryBase):
    """
    xTB theory factory.

    Config:
        {
            "method": "GFN2",        # xTB method (GFN1, GFN2, GFN0)
            "runmode": "inputfile"   # Run mode (inputfile, etc.)
        }
    """

    @property
    def name(self) -> str:
        return "xtb"

    def validate_config(self) -> bool:
        """Validate xTB configuration."""
        required = ["method", "runmode"]
        return all(k in self.config for k in required)

    def create_theory(
        self, numcores: int = 1, printlevel: int = 1, **kwargs
    ) -> Any:
        """
        Create xTB theory instance.

        Args:
            numcores: Number of CPU cores
            printlevel: Print level for output

        Returns:
            xTBTheory instance
        """
        from ash import xTBTheory

        if not self.validate_config():
            raise ValueError(f"Invalid xTB config: {self.config}")

        return xTBTheory(
            xtbmethod=self.config["method"],
            runmode=self.config["runmode"],
            numcores=numcores,
            printlevel=printlevel,
            filename="xtb_",
        )

    def get_default_config(self) -> Dict[str, Any]:
        """Return default xTB configuration."""
        return {
            "method": "GFN2",
            "runmode": "inputfile",
        }
