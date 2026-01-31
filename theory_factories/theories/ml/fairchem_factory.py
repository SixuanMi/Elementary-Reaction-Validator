#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FairChem theory factory.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from ...base import TheoryFactoryBase
from ...registry import register_theory


@register_theory("fairchem")
class FairChemFactory(TheoryFactoryBase):
    """
    FairChem theory factory.

    Config:
        {
            "model_file": "path/to/uma-s-1p1.pt",  # Path to FairChem model file
            "task_name": "omol"                    # Task name (oc20, omol, omat, odac, omc)
            "device": "cpu"                        # Device (cpu, cuda, opencl, mps)
        }
    """

    @property
    def name(self) -> str:
        return "fairchem"

    def validate_config(self) -> bool:
        """Validate FairChem configuration."""
        if "model_file" not in self.config:
            return False
        if "task_name" not in self.config:
            return False
        model_file = Path(self.config["model_file"])
        return model_file.exists()

    def create_theory(
        self, numcores: int = 1, printlevel: int = 1, **kwargs
    ) -> Any:
        """
        Create FairChem theory instance.

        Args:
            numcores: Number of CPU cores
            printlevel: Print level for output

        Returns:
            FairchemTheory instance
        """
        from ash import FairchemTheory

        if not self.validate_config():
            raise ValueError(f"Invalid FairChem config: {self.config}")

        return FairchemTheory(
            model_file=self.config["model_file"],
            task_name=self.config["task_name"],
            device=self.config.get("device", "cpu"),
            numcores=numcores,
        )

    def get_default_config(self) -> Dict[str, Any]:
        """Return default FairChem configuration."""
        return {
            "model_file": "",
            "task_name": "omol",
            "device": "cpu",
        }
