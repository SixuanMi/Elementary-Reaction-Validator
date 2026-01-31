#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MACE theory factory.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from ...base import TheoryFactoryBase
from ...registry import register_theory


@register_theory("mace")
class MACEFactory(TheoryFactoryBase):
    """
    MACE theory factory.

    Config:
        {
            "model_file": "path/to/model.model",  # Path to MACE model file
            "device": "cpu"                        # Device (cpu, cuda, etc.)
        }
    """

    @property
    def name(self) -> str:
        return "mace"

    def validate_config(self) -> bool:
        """Validate MACE configuration."""
        if "model_file" not in self.config:
            return False
        model_file = Path(self.config["model_file"])
        return model_file.exists()

    def create_theory(
        self, numcores: int = 1, printlevel: int = 1, **kwargs
    ) -> Any:
        """
        Create MACE theory instance.

        Args:
            numcores: Number of CPU cores
            printlevel: Print level for output

        Returns:
            MACETheory instance
        """
        from ash import MACETheory

        if not self.validate_config():
            raise ValueError(f"Invalid MACE config: {self.config}")

        return MACETheory(
            model_file=self.config["model_file"],
            device=self.config.get("device", "cpu"),
            numcores=numcores,
            printlevel=printlevel,
        )

    def get_default_config(self) -> Dict[str, Any]:
        """Return default MACE configuration."""
        return {
            "model_file": "",
            "device": "cpu",
        }
