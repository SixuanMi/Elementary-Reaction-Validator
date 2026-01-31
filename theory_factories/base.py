#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base classes for theory factories.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional


class TheoryFactoryBase(ABC):
    """
    Abstract base class for theory factory implementations.

    Subclasses must implement:
    - name: Return the theory name
    - create_theory(): Create a theory instance
    - validate_config(): Validate the configuration
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the factory with optional configuration.

        Args:
            config: Method-specific configuration parameters
        """
        self.config = config or {}

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the theory name."""
        pass

    @abstractmethod
    def create_theory(self, **kwargs) -> Any:
        """
        Create a theory instance.

        Common parameters (passed via kwargs):
        - numcores: Number of CPU cores
        - printlevel: Print level for theory output

        Returns:
            Theory instance (xTBTheory, ORCATheory, MACETheory, etc.)
        """
        pass

    @abstractmethod
    def validate_config(self) -> bool:
        """
        Validate the configuration.

        Returns:
            True if configuration is valid, False otherwise
        """
        pass

    def get_default_config(self) -> Dict[str, Any]:
        """
        Return the default configuration for this theory.

        Returns:
            Dictionary with default configuration values
        """
        return {}
