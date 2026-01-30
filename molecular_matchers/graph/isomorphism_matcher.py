#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Graph isomorphism-based molecular matcher using Pymatgen + NetworkX.

Uses molecular graph representation with bond connectivity to determine
if two structures are topologically equivalent (isomorphic).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from ..base import MoleculeMatcher
from ..registry import register_matcher


@register_matcher("graph_isomorphism")
class GraphIsomorphismMatcher(MoleculeMatcher):
    """
    Graph isomorphism matcher using Pymatgen MoleculeGraph + NetworkX.

    Config:
        local_env: str = "OpenBabelNN"  # Local environment strategy for bond determination

    Uses OpenBabelNN (OpenBabel neighbor finding) to determine molecular graph,
    then NetworkX isomorphism algorithms for comparison.
    """

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self._local_env = self.config.get("local_env", "OpenBabelNN")
        self._mol_graphs = {}  # Cache for MoleculeGraph objects

    @property
    def name(self) -> str:
        return "graph_isomorphism"

    def _check_dependencies(self) -> bool:
        """Check if Pymatgen and NetworkX are available."""
        try:
            from networkx.algorithms import isomorphism as iso
            from pymatgen.core import Molecule
            from pymatgen.analysis.graphs import MoleculeGraph
            from pymatgen.analysis.local_env import OpenBabelNN
            return True
        except Exception:
            return False

    def _get_local_env_strategy(self):
        """Get the local environment strategy class."""
        from pymatgen.analysis.local_env import OpenBabelNN, CovalentBondNN, MinimumDistanceNN

        strategies = {
            "OpenBabelNN": OpenBabelNN,
            "CovalentBondNN": CovalentBondNN,
            "MinimumDistanceNN": MinimumDistanceNN,
        }
        return strategies.get(self._local_env, OpenBabelNN)()

    def _create_molecule_graph(self, xyz_file: Path):
        """
        Create a Pymatgen MoleculeGraph from XYZ file.

        Args:
            xyz_file: Path to XYZ structure file

        Returns:
            MoleculeGraph object, or None if creation fails
        """
        if not self._check_dependencies():
            return None

        try:
            from pymatgen.core import Molecule
            from pymatgen.analysis.graphs import MoleculeGraph

            # Use cached graph if available
            cache_key = str(xyz_file)
            if cache_key in self._mol_graphs:
                return self._mol_graphs[cache_key]

            # Create Molecule from XYZ
            mol = Molecule.from_file(str(xyz_file))

            # Create MoleculeGraph with local environment strategy
            # Use new API: from_local_env_strategy (old with_local_env_strategy is deprecated)
            local_env_strategy = self._get_local_env_strategy()
            try:
                mol_graph = MoleculeGraph.from_local_env_strategy(mol, local_env_strategy)
            except AttributeError:
                # Fallback to old API for older pymatgen versions
                mol_graph = MoleculeGraph.with_local_env_strategy(mol, local_env_strategy)

            # Cache for reuse
            self._mol_graphs[cache_key] = mol_graph
            return mol_graph

        except Exception:
            return None

    def compute_signature(self, xyz_file: Path) -> Optional[Any]:
        """
        Create molecular graph signature from XYZ file.

        Args:
            xyz_file: Path to XYZ structure file

        Returns:
            MoleculeGraph object for isomorphism comparison, or None if creation fails
        """
        return self._create_molecule_graph(xyz_file)

    def compare(self, sig1: Any, sig2: Any) -> bool:
        """
        Compare two molecular graphs for isomorphism.

        Args:
            sig1: First MoleculeGraph
            sig2: Second MoleculeGraph

        Returns:
            True if graphs are isomorphic (same connectivity)
        """
        if sig1 is None or sig2 is None:
            return False

        try:
            from networkx.algorithms import isomorphism as iso

            gm = iso.GraphMatcher(sig1.graph, sig2.graph)
            return gm.is_isomorphic()
        except Exception:
            return False

    def cross_isomorphism(
        self,
        true_r: Path,
        true_p: Path,
        irc_r: Path,
        irc_p: Path,
    ) -> Any:
        """
        Cross-compare IRC endpoints with original reactant/product using graph isomorphism.

        This override adds the original pymatgen-specific output fields
        (reactant_isomorphic, product_isomorphic, both_isomorphic) for compatibility.

        Args:
            true_r: Original reactant XYZ file
            true_p: Original product XYZ file
            irc_r: IRC backward endpoint XYZ file
            irc_p: IRC forward endpoint XYZ file

        Returns:
            MatchResult with isomorphism results
        """
        from ..base import MatchResult

        try:
            # Create molecular graphs
            mg_true_r = self._create_molecule_graph(true_r)
            mg_true_p = self._create_molecule_graph(true_p)
            mg_irc_r = self._create_molecule_graph(irc_r)
            mg_irc_p = self._create_molecule_graph(irc_p)

            # Check for failures
            if None in (mg_true_r, mg_true_p, mg_irc_r, mg_irc_p):
                return MatchResult(
                    performed=True,
                    success=False,
                    endpoint_match="Error",
                    method=self.name,
                    error="Failed to create one or more molecular graphs",
                )

            # NetworkX isomorphism checks - all 4 combinations
            from networkx.algorithms import isomorphism as iso

            gm_tr = iso.GraphMatcher(mg_true_r.graph, mg_irc_r.graph)
            gm_tp = iso.GraphMatcher(mg_true_p.graph, mg_irc_p.graph)
            gm_cross_rp = iso.GraphMatcher(mg_true_r.graph, mg_irc_p.graph)
            gm_cross_pr = iso.GraphMatcher(mg_true_p.graph, mg_irc_r.graph)

            direct_match_r = gm_tr.is_isomorphic()
            direct_match_p = gm_tp.is_isomorphic()
            cross_match_rp = gm_cross_rp.is_isomorphic()
            cross_match_pr = gm_cross_pr.is_isomorphic()

            reactant_iso = direct_match_r
            product_iso = direct_match_p
            both_iso = direct_match_r and direct_match_p

            # Build matches tuple (true_r==irc_r, true_p==irc_p, true_r==irc_p, true_p==irc_r)
            matches = (direct_match_r, direct_match_p, cross_match_rp, cross_match_pr)

            # Determine endpoint match category
            if matches in [(True, True, False, False), (False, False, True, True)]:
                endpoint_match = "2-end match"
            elif any(matches):
                endpoint_match = "1-end match"
            else:
                endpoint_match = "No match"

            # Determine reaction type (chemical vs conformational)
            try:
                gm_irc = iso.GraphMatcher(mg_irc_r.graph, mg_irc_p.graph)
                rxn_status = "Conformational change" if gm_irc.is_isomorphic() else "Chemical reaction"
            except Exception:
                rxn_status = "Chemical reaction"

            return MatchResult(
                performed=True,
                success=True,
                endpoint_match=endpoint_match,
                matches=matches,
                rxn_status=rxn_status,
                method=self.name,
                metadata={
                    "reactant_isomorphic": direct_match_r,
                    "product_isomorphic": direct_match_p,
                    "both_isomorphic": both_iso,
                    "cross_match_rp": cross_match_rp,
                    "cross_match_pr": cross_match_pr,
                    "local_env_strategy": self._local_env,
                },
            )

        except Exception as e:
            return MatchResult(
                performed=True,
                success=False,
                endpoint_match="Error",
                method=self.name,
                error=str(e),
            )
