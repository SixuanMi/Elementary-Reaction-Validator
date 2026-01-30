#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matcher registry for dynamic method selection.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Type

from .base import MoleculeMatcher


class MatcherRegistry:
    """
    Global registry for molecule matcher implementations.
    """

    _matchers: Dict[str, Type[MoleculeMatcher]] = {}

    @classmethod
    def register(cls, name: str, matcher_class: Type[MoleculeMatcher]) -> None:
        """Register a matcher class under a name."""
        cls._matchers[name] = matcher_class

    @classmethod
    def get(cls, name: str) -> Optional[Type[MoleculeMatcher]]:
        """Get a matcher class by name."""
        return cls._matchers.get(name)

    @classmethod
    def list_names(cls) -> List[str]:
        """List all registered matcher names."""
        return sorted(cls._matchers.keys())

    @classmethod
    def create(
        cls, name: str, config: Optional[Dict[str, Any]] = None
    ) -> Optional[MoleculeMatcher]:
        """Create a matcher instance by name with optional config."""
        matcher_class = cls.get(name)
        if matcher_class is None:
            return None
        return matcher_class(config=config)


def register_matcher(name: str):
    """
    Decorator to register a matcher class.

    Usage:
        @register_matcher("smiles_openbabel")
        class OpenBabelSmilesMatcher(MoleculeMatcher):
            ...
    """

    def decorator(matcher_class: Type[MoleculeMatcher]):
        MatcherRegistry.register(name, matcher_class)
        return matcher_class

    return decorator


def get_matcher(name: str, config: Optional[Dict[str, Any]] = None) -> Optional[MoleculeMatcher]:
    """Create a matcher instance by name."""
    return MatcherRegistry.create(name, config=config)


def list_matchers() -> List[str]:
    """List all registered matcher names."""
    return MatcherRegistry.list_names()


def create_matcher_from_config(config: Dict[str, Any]) -> Optional[MoleculeMatcher]:
    """
    Create a matcher instance from configuration dict.

    Config format:
        {
            "method": "smiles_openbabel",
            "smiles_openbabel": {...},  # method-specific config
            ...
        }

    The method-specific config section name must match the method name.
    """
    method = config.get("method")
    if not method:
        return None

    # Extract method-specific config
    method_config = config.get(method, {})

    return get_matcher(method, config=method_config)
