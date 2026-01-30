#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Molecular Matchers Module

Provides a pluggable architecture for molecular structure comparison/matching.
Supports multiple matching strategies:
- SMILES-based (OpenBabel, RDKit)
- Graph isomorphism (Pymatgen + NetworkX)
- 3D geometry (RMSD)
- Descriptor similarity (SOAP)
"""

from .base import MoleculeMatcher, MatchResult
from .registry import (
    get_matcher,
    list_matchers,
    register_matcher,
    create_matcher_from_config,
)

# Import all matchers to register them
from .smiles import OpenBabelSmilesMatcher, RDKitSmilesMatcher
from .graph import GraphIsomorphismMatcher
from .geometry import RMSDMatcher
from .descriptor import SOAPMatcher

__all__ = [
    # Base classes
    "MoleculeMatcher",
    "MatchResult",
    # Registry functions
    "get_matcher",
    "list_matchers",
    "register_matcher",
    "create_matcher_from_config",
    # Concrete matchers
    "OpenBabelSmilesMatcher",
    "RDKitSmilesMatcher",
    "GraphIsomorphismMatcher",
    "RMSDMatcher",
    "SOAPMatcher",
]
