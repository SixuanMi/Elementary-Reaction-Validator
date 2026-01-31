#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Theory factory registry for dynamic theory selection.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Type

from .base import TheoryFactoryBase


class TheoryRegistry:
    """
    Global registry for theory factory implementations.
    """

    _factories: Dict[str, Type[TheoryFactoryBase]] = {}

    @classmethod
    def register(cls, name: str, factory_class: Type[TheoryFactoryBase]) -> None:
        """Register a factory class under a name."""
        cls._factories[name] = factory_class

    @classmethod
    def get(cls, name: str) -> Optional[Type[TheoryFactoryBase]]:
        """Get a factory class by name."""
        return cls._factories.get(name)

    @classmethod
    def list_names(cls) -> List[str]:
        """List all registered factory names."""
        return sorted(cls._factories.keys())

    @classmethod
    def create(
        cls, name: str, config: Optional[Dict[str, Any]] = None
    ) -> Optional[TheoryFactoryBase]:
        """Create a factory instance by name with optional config."""
        factory_class = cls.get(name)
        if factory_class is None:
            return None
        return factory_class(config=config)


def register_theory(name: str):
    """
    Decorator to register a theory factory class.

    Usage:
        @register_theory("xtb")
        class XTBFactory(TheoryFactoryBase):
            ...
    """

    def decorator(factory_class: Type[TheoryFactoryBase]):
        TheoryRegistry.register(name, factory_class)
        return factory_class

    return decorator


def get_theory_factory(name: str, config: Optional[Dict[str, Any]] = None) -> Optional[TheoryFactoryBase]:
    """Create a factory instance by name."""
    return TheoryRegistry.create(name, config=config)


def list_theories() -> List[str]:
    """List all registered theory names."""
    return TheoryRegistry.list_names()


def create_theory_from_config(config: Dict[str, Any]) -> Optional[Any]:
    """
    Create a theory instance from configuration dict.

    Config format:
        {
            "type": "mace",           # theory type (required)
            "numcores": 1,             # common params (optional)
            "printlevel": 1,           # common params (optional)
            "mace": {                  # method-specific config (optional)
                "model_file": "path/to/model.model",
                "device": "cpu"
            }
        }

    The method-specific config section name must match the theory type.
    For example, if type="mace", the method config should be under "mace" key.

    Args:
        config: Configuration dictionary

    Returns:
        Theory instance (xTBTheory, ORCATheory, MACETheory, etc.)
        Returns None if theory type is not found
    """
    theory_type = config.get("type")
    if not theory_type:
        return None

    # Extract common parameters
    common_params = {
        "numcores": config.get("numcores", 1),
        "printlevel": config.get("printlevel", 1),
    }

    # Extract method-specific config (key must match theory type)
    method_config = config.get(theory_type, {})

    factory = get_theory_factory(theory_type, config=method_config)
    if factory is None:
        return None

    return factory.create_theory(**common_params)
