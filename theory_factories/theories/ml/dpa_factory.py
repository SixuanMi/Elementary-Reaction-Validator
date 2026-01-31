#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DPA3 theory factory.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from ...base import TheoryFactoryBase
from ...registry import register_theory


@register_theory("dpa")
class DPAFactory(TheoryFactoryBase):
    """
    DPA3 (Deep Potential Approximation) theory factory.

    Config:
        {
            "model_file": "path/to/DPA-3.1-3M.pt",  # Path to DPA model file
            "head": "Organic_Reactions",             # Model head for multitask models
            "device": "cpu"                          # (placeholder, for compatibility)
        }
    """

    @property
    def name(self) -> str:
        return "dpa"

    def validate_config(self) -> bool:
        """Validate DPA configuration."""
        if "model_file" not in self.config:
            return False
        model_file = Path(self.config["model_file"])
        return model_file.exists()

    def create_theory(
        self, numcores: int = 1, printlevel: int = 1, **kwargs
    ) -> Any:
        """
        Create DPA3 theory instance.

        Args:
            numcores: Number of CPU cores
            printlevel: Print level for output

        Returns:
            DPA3Theory instance
        """
        # Import from local interface file
        import sys
        repo_root = Path(__file__).parent.parent.parent.parent.parent
        sys.path.insert(0, str(repo_root))
        from interface_DPA3 import DPA3Theory

        if not self.validate_config():
            raise ValueError(f"Invalid DPA config: {self.config}")

        return DPA3Theory(
            DPAdir=self.config["model_file"],
            head=self.config.get("head"),
            numcores=numcores,
        )

    def get_default_config(self) -> Dict[str, Any]:
        """Return default DPA configuration."""
        return {
            "model_file": "",
            "head": "Organic_Reactions",  # Default head for organic reactions
            "device": "cpu",  # Placeholder for compatibility
        }
