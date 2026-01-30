#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base classes for molecular matchers.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class MatchResult:
    """
    Result of a molecular matching operation.

    Attributes:
        performed: Whether the matching was attempted
        success: Whether the matching succeeded without errors
        endpoint_match: Match category ("2-end match", "1-end match", "No match", "Error")
        matches: Tuple of 4 booleans (true_r==irc_r, true_p==irc_p, true_r==irc_p, true_p==irc_r)
        rxn_status: "Chemical reaction" or "Conformational change"
        method: Name of the matching method used
        metadata: Additional method-specific data
        error: Error message if success=False
    """

    performed: bool = True
    success: bool = True
    endpoint_match: str = "No match"
    matches: Tuple[bool, bool, bool, bool] = (False, False, False, False)
    rxn_status: str = "Chemical reaction"
    method: str = "unknown"
    metadata: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "performed": self.performed,
            "success": self.success,
            "endpoint_match": self.endpoint_match,
            "matches": self.matches,
            "rxn_status": self.rxn_status,
            "method": self.method,
            "metadata": self.metadata,
            "error": self.error,
        }


class MoleculeMatcher(ABC):
    """
    Abstract base class for molecular structure matching.

    Subclasses must implement:
    - compute_signature(): Compute molecular signature/descriptor from XYZ
    - compare(): Compare two signatures for equality/similarity
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the matcher with optional configuration.

        Args:
            config: Method-specific configuration parameters
        """
        self.config = config or {}
        self._name: Optional[str] = None

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of this matcher."""
        pass

    @abstractmethod
    def compute_signature(self, xyz_file: Path) -> Any:
        """
        Compute molecular signature/descriptor from XYZ file.

        Args:
            xyz_file: Path to XYZ structure file

        Returns:
            Molecular signature (type depends on matcher)
            Returns None if computation fails

        Raises:
            Exception: If computation fails (will be caught by caller)
        """
        pass

    @abstractmethod
    def compare(self, sig1: Any, sig2: Any) -> bool:
        """
        Compare two molecular signatures for equality/similarity.

        Args:
            sig1: First signature (from compute_signature)
            sig2: Second signature (from compute_signature)

        Returns:
            True if signatures match, False otherwise
        """
        pass

    def are_same_molecule(self, sig1: Any, sig2: Any) -> bool:
        """
        Check if two signatures represent the same molecule.
        Default implementation uses compare().

        Override this if you need different logic for "same molecule"
        vs "endpoint match" (e.g., for conformational changes).
        """
        return self.compare(sig1, sig2)

    def cross_isomorphism(
        self,
        true_r: Path,
        true_p: Path,
        irc_r: Path,
        irc_p: Path,
    ) -> MatchResult:
        """
        Cross-compare IRC endpoints with original reactant/product.

        This is the main entry point for endpoint matching validation.

        Args:
            true_r: Original reactant XYZ file
            true_p: Original product XYZ file
            irc_r: IRC backward endpoint XYZ file
            irc_p: IRC forward endpoint XYZ file

        Returns:
            MatchResult with comparison results
        """
        try:
            # Compute signatures for all 4 structures
            sig_true_r = self.compute_signature(true_r)
            sig_true_p = self.compute_signature(true_p)
            sig_irc_r = self.compute_signature(irc_r)
            sig_irc_p = self.compute_signature(irc_p)

            # Check for computation failures
            if None in (sig_true_r, sig_true_p, sig_irc_r, sig_irc_p):
                return MatchResult(
                    performed=True,
                    success=False,
                    endpoint_match="Error",
                    method=self.name,
                    metadata={
                        "error_detail": "Failed to compute one or more signatures",
                    },
                    error="Failed to compute one or more signatures",
                )

            # Cross-compare
            matches = (
                self.compare(sig_true_r, sig_irc_r),
                self.compare(sig_true_p, sig_irc_p),
                self.compare(sig_true_r, sig_irc_p),
                self.compare(sig_true_p, sig_irc_r),
            )

            # Determine reaction type (chemical vs conformational)
            rxn_status = (
                "Conformational change"
                if self.are_same_molecule(sig_irc_r, sig_irc_p)
                else "Chemical reaction"
            )

            # Determine endpoint match category
            if matches in [(True, True, False, False), (False, False, True, True)]:
                endpoint_match = "2-end match"
            elif any(matches):
                endpoint_match = "1-end match"
            else:
                endpoint_match = "No match"

            return MatchResult(
                performed=True,
                success=True,
                endpoint_match=endpoint_match,
                matches=matches,
                rxn_status=rxn_status,
                method=self.name,
                metadata=self._build_cross_isomorphism_metadata(
                    sig_true_r, sig_true_p, sig_irc_r, sig_irc_p
                ),
            )

        except Exception as e:
            return MatchResult(
                performed=True,
                success=False,
                endpoint_match="Error",
                method=self.name,
                error=str(e),
            )

    def _serialize_signature(self, sig: Any) -> Any:
        """Serialize signature for JSON output."""
        return sig

    def _build_cross_isomorphism_metadata(
        self, sig_true_r: Any, sig_true_p: Any, sig_irc_r: Any, sig_irc_p: Any
    ) -> Dict[str, Any]:
        """
        Build metadata dict for cross-isomorphism results.
        Override in subclasses to add method-specific metadata.
        """
        return {
            "sig_true_r": self._serialize_signature(sig_true_r),
            "sig_true_p": self._serialize_signature(sig_true_p),
            "sig_irc_r": self._serialize_signature(sig_irc_r),
            "sig_irc_p": self._serialize_signature(sig_irc_p),
        }
