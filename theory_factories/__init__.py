#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Theory Factories Module

Provides a pluggable architecture for theory object creation.
Supports multiple theory types:
- Semi-empirical: xTB
- DFT: ORCA
- Machine Learning: MACE, FairChem
"""

from __future__ import annotations

from .base import TheoryFactoryBase
from .registry import (
    TheoryRegistry,
    get_theory_factory,
    list_theories,
    register_theory,
    create_theory_from_config,
)

# Import all factories to register them
from .theories import XTBFactory, ORCAFactory, MACEFactory, FairChemFactory

__all__ = [
    # Base classes
    "TheoryFactoryBase",
    # Registry classes/functions
    "TheoryRegistry",
    "get_theory_factory",
    "list_theories",
    "register_theory",
    "create_theory_from_config",
    # Concrete factories
    "XTBFactory",
    "ORCAFactory",
    "MACEFactory",
    "FairChemFactory",
]
