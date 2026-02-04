#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import contextlib
import json
import os
import pickle
import re
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml

# Suppress matplotlib fontManager warning that can interfere with geometric log parsing
# This warning appears as "generated new fontManager" in logs and can break step parsing
import warnings
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.font_manager
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


def _add_local_ash_to_syspath() -> Tuple[Path, Path]:
    repo_root = Path(__file__).resolve().parent
    local_ash_root = repo_root / "ash"
    sys.path.insert(0, str(local_ash_root))
    # Also add repo_root to sys.path for local modules like interfaces
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    return repo_root, local_ash_root


@contextlib.contextmanager
def _pushd(path: Path):
    old_cwd = Path.cwd()
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_cwd)


def _list_files(root: Path) -> List[str]:
    return sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file())


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    return True


@dataclass
class StageSummary:
    name: str
    workdir: str
    ok: bool
    result: Dict[str, Any]
    file_count: int
    files_txt: str
    error: Optional[str] = None


def _run_stage_inventory(stage_dir: Path, stage_name: str) -> Tuple[int, Path]:
    files = _list_files(stage_dir)
    inventory_path = stage_dir / f"{stage_name}_files.txt"
    _write_text(inventory_path, "\n".join(files) + ("\n" if files else ""))
    return len(files), inventory_path


def _as_float_or_none(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except Exception:
        return None


_ANSI_RE = re.compile(r"\x1B\[[0-?]*[ -/]*[@-~]")


def _strip_ansi(text: str) -> str:
    return _ANSI_RE.sub("", text)


def _read_json_if_exists(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _file_contains_text(path: Path, needle: str) -> bool:
    if not path.exists():
        return False
    try:
        return needle in path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return False


def _geometric_converged(log_path: Path) -> bool:
    # geomeTRIC prints this exact string on convergence.
    return _file_contains_text(log_path, "Converged! =D")


def _count_significant_imaginary(freqs: List[float], imag_threshold_cm1: float) -> Tuple[int, List[float]]:
    imag = [float(f) for f in freqs if float(f) < -float(imag_threshold_cm1)]
    return len(imag), imag


def _serialize_mep_energies(mep: Any) -> Optional[Dict[str, Optional[float]]]:
    if not isinstance(mep, dict):
        return None
    out: Dict[str, Optional[float]] = {}
    for k, v in mep.items():
        out[str(k)] = _as_float_or_none(v)
    return out




_GEOM_STEP_RE = re.compile(r"^Step\s+(\d+)\s*:")
_XYZ_ITER_RE = re.compile(r"\bIteration\s+(\d+)\b")
_XYZ_ENERGY_RE = re.compile(r"\bEnergy\s+(-?\d+(?:\.\d+)?)\b")


def _parse_geometric_irc_end_steps(log_path: Path) -> Tuple[Optional[int], Optional[int], Dict[str, Any]]:
    """
    Parse geomeTRIC IRC log to find the final step index for:
      - forward direction (ending right before 'IRC forward direction converged')
      - backward direction (ending at the last 'Step N :' after 'IRC backward direction starts here')
    """
    meta: Dict[str, Any] = {
        "has_forward_converged_line": False,
        "has_backward_start_line": False,
        "has_final_converged_line": False,
    }
    if not log_path.exists():
        return None, None, meta

    text = _strip_ansi(log_path.read_text(encoding="utf-8", errors="ignore"))
    lines = text.splitlines()

    idx_forward_conv = None
    idx_backward_start = None
    for i, ln in enumerate(lines):
        if "IRC forward direction converged" in ln:
            idx_forward_conv = i
            meta["has_forward_converged_line"] = True
        if "IRC backward direction starts here" in ln:
            idx_backward_start = i
            meta["has_backward_start_line"] = True
    meta["has_final_converged_line"] = "Converged! =D" in text

    forward_end: Optional[int] = None
    backward_end: Optional[int] = None

    if idx_forward_conv is not None:
        for ln in lines[:idx_forward_conv]:
            m = _GEOM_STEP_RE.match(ln.strip())
            if m:
                forward_end = int(m.group(1))

    if idx_backward_start is not None:
        for ln in lines[idx_backward_start:]:
            m = _GEOM_STEP_RE.match(ln.strip())
            if m:
                backward_end = int(m.group(1))

    return forward_end, backward_end, meta




def _parse_xyz_blocks_in_order(xyz_path: Path) -> List[List[str]]:
    """
    Parse an XYZ trajectory into a list of blocks, preserving file order.
    Each block is [natoms_line, comment_line, atom_lines...].
    """
    text = xyz_path.read_text(encoding="utf-8", errors="ignore")
    lines = [ln.rstrip("\n") for ln in text.splitlines()]

    blocks: List[List[str]] = []
    idx = 0
    while idx < len(lines):
        if not lines[idx].strip():
            idx += 1
            continue
        natoms = int(lines[idx].split()[0])
        if idx + 1 + natoms >= len(lines):
            raise ValueError(f"Truncated XYZ block starting at line {idx + 1} in {xyz_path}")
        blocks.append(lines[idx : idx + 2 + natoms])
        idx += 2 + natoms
    return blocks


def _extract_geometric_irc_endpoints(
    irc_log_path: Path,
    irc_traj_xyz: Path,
    out_forward_xyz: Path,
    out_backward_xyz: Path,
) -> Dict[str, Any]:
    """
    Extract IRC forward/backward endpoint geometries from geomeTRIC outputs.
    """
    forward_end, backward_end, meta = _parse_geometric_irc_end_steps(irc_log_path)
    if not irc_traj_xyz.exists():
        raise FileNotFoundError(f"Missing IRC trajectory XYZ: {irc_traj_xyz}")

    # NOTE: geomeTRIC IRC trajectories can contain duplicated "Iteration N" labels
    # (e.g. 48..0 then 0..72), so indexing by iteration loses information.
    # For robust endpoint extraction, use the first and last XYZ frames.
    blocks = _parse_xyz_blocks_in_order(irc_traj_xyz)
    if len(blocks) < 2:
        raise ValueError(f"IRC trajectory has <2 frames: {irc_traj_xyz}")

    first_block = blocks[0]
    last_block = blocks[-1]

    def _iter_from_block(block: List[str]) -> Optional[int]:
        if len(block) < 2:
            return None
        m = _XYZ_ITER_RE.search(block[1])
        return int(m.group(1)) if m else None

    def _energy_from_block(block: List[str]) -> Optional[float]:
        if len(block) < 2:
            return None
        m = _XYZ_ENERGY_RE.search(block[1])
        return float(m.group(1)) if m else None

    first_it = _iter_from_block(first_block)
    last_it = _iter_from_block(last_block)

    forward_block = first_block
    backward_block = last_block
    selection = "first_last"
    if forward_end is not None and backward_end is not None and first_it is not None and last_it is not None:
        if first_it == forward_end and last_it == backward_end:
            selection = "first_last_matched_log_steps"
        elif first_it == backward_end and last_it == forward_end:
            forward_block = last_block
            backward_block = first_block
            selection = "first_last_swapped_to_match_log_steps"
        else:
            selection = "first_last_log_steps_mismatch"

    out_forward_xyz.parent.mkdir(parents=True, exist_ok=True)
    out_forward_xyz.write_text("\n".join(forward_block) + "\n", encoding="utf-8")
    out_backward_xyz.parent.mkdir(parents=True, exist_ok=True)
    out_backward_xyz.write_text("\n".join(backward_block) + "\n", encoding="utf-8")

    return {
        "selection_strategy": selection,
        "first_frame_iteration": first_it,
        "last_frame_iteration": last_it,
        "log_forward_end_step": forward_end,
        "log_backward_end_step": backward_end,
        "forward_end_step": _iter_from_block(forward_block),
        "backward_end_step": _iter_from_block(backward_block),
        "forward_end_energy": _energy_from_block(forward_block),
        "backward_end_energy": _energy_from_block(backward_block),
        "log_meta": meta,
        "out_forward_xyz": str(out_forward_xyz),
        "out_backward_xyz": str(out_backward_xyz),
    }


def _write_workflow_summary(run_dir: Path, summary: Dict[str, Any], stages: List[StageSummary]) -> Path:
    summary["stages"] = [s.__dict__ for s in stages]
    out_path = run_dir / "workflow_summary.json"
    _write_text(out_path, json.dumps(summary, indent=2, ensure_ascii=False))
    return out_path


# ============================================================================
# Molecular Matcher Module (cross-isomorphism with multiple methods)
# ============================================================================

# Import the modular molecular matchers
# This provides pluggable endpoint matching strategies:
# - smiles_openbabel: OpenBabel SMILES comparison
# - smiles_rdkit: RDKit SMILES comparison with multiple bond strategies
# - graph_isomorphism: Pymatgen + NetworkX graph isomorphism
# - rmsd: 3D RMSD structural comparison
# - soap: SOAP descriptor similarity
def _get_endpoint_matcher(endpoint_match_config: Dict[str, Any]):
    """
    Create a molecular matcher instance from config.

    Args:
        endpoint_match_config: Config dict with 'method' and method-specific settings

    Returns:
        MoleculeMatcher instance or None if method not found
    """
    from molecular_matchers import create_matcher_from_config

    return create_matcher_from_config(endpoint_match_config)


def _cross_isomorphism(
    true_r: Path,
    true_p: Path,
    irc_r: Path,
    irc_p: Path,
    matcher=None,
) -> Dict[str, Any]:
    """
    Cross-compare endpoints using the configured molecular matcher.

    This is a backward-compatible wrapper that uses the new modular matcher system.
    If no matcher is provided, falls back to the default OpenBabel SMILES method.

    Args:
        true_r: Original reactant XYZ file
        true_p: Original product XYZ file
        irc_r: IRC backward endpoint XYZ file
        irc_p: IRC forward endpoint XYZ file
        matcher: Optional MoleculeMatcher instance (if None, uses OpenBabel default)

    Returns:
        Dict with match results
    """
    # Use default OpenBabel matcher if none provided
    if matcher is None:
        from molecular_matchers import OpenBabelSmilesMatcher
        matcher = OpenBabelSmilesMatcher()

    # Run cross-isomorphism check
    result = matcher.cross_isomorphism(true_r, true_p, irc_r, irc_p)
    result_dict = result.to_dict()

    # Flatten metadata for backward compatibility (only for SMILES methods)
    # For other methods, keep metadata structured to avoid duplication
    metadata = result_dict.get("metadata", {})

    if "smiles_openbabel" in result.method or "smiles_rdkit" in result.method:
        # For SMILES, flatten smi_* fields to top level for convenience
        result_dict["smi_true_r"] = metadata.get("sig_true_r", "")
        result_dict["smi_true_p"] = metadata.get("sig_true_p", "")
        result_dict["smi_irc_r"] = metadata.get("sig_irc_r", "")
        result_dict["smi_irc_p"] = metadata.get("sig_irc_p", "")
        # Remove from metadata to avoid duplication
        for key in ["sig_true_r", "sig_true_p", "sig_irc_r", "sig_irc_p"]:
            metadata.pop(key, None)

    # For other methods, keep metadata structured (no flattening)
    # This keeps the JSON clean and avoids duplication

    return result_dict


# Legacy function kept for backward compatibility
def _canonical_smiles_from_xyz(xyz_file: Path) -> str:
    """
    Legacy function: Best-effort XYZ -> canonical SMILES using OpenBabel/pybel.
    Returns "" on failure.

    Note: This is kept for backward compatibility.
    New code should use the molecular_matchers module directly.
    """
    from molecular_matchers import OpenBabelSmilesMatcher

    matcher = OpenBabelSmilesMatcher()
    result = matcher.compute_signature(xyz_file)
    return result if result else ""


def _read_xyz_elements_coords(xyz_path: Path) -> Tuple[List[str], np.ndarray]:
    """Read XYZ file, return (elements, coords_array_in_Angstrom)"""
    with open(xyz_path) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    elements = []
    coords = []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        elements.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elements, np.array(coords)


def _read_energy_file(energy_file: Path) -> Optional[float]:
    """Read ASH energy file, return energy in Eh"""
    if not energy_file.exists():
        return None
    match = re.search(r"-?\d+\.\d+", energy_file.read_text())
    return float(match.group()) if match else None


def _read_gradient_file(grad_file: Path) -> Optional[np.ndarray]:
    """Read ASH gradient file, return (N,3) array in Eh/Angstrom"""
    if not grad_file.exists():
        return None
    grad_data = []
    for line in grad_file.read_text().splitlines():
        if "$grad" in line or "$end" in line or "cycle" in line:
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                grad_data.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                pass
    return np.array(grad_data) if grad_data else None


def _read_irc_trajectory(irc_traj_xyz: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Parse IRC trajectory XYZ, return (coords_array, energies_array)

    Format:
        12
        Iteration 48 Energy -21.56784433
        C  x  y  z
        ...
    """
    if not irc_traj_xyz.exists():
        return None, None

    with open(irc_traj_xyz) as f:
        lines = f.readlines()

    frames = []
    energies = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            natoms = int(lines[i].strip())
        except ValueError:
            i += 1
            continue

        if i + 1 >= len(lines):
            break
        comment = lines[i + 1].strip()

        energy_match = re.search(r"Energy\s+(-?\d+\.\d+)", comment)
        energy = float(energy_match.group(1)) if energy_match else None

        coords = []
        for j in range(i + 2, min(i + 2 + natoms, len(lines))):
            parts = lines[j].split()
            if len(parts) >= 4:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        if len(coords) == natoms:
            frames.append(np.array(coords))
            energies.append(energy)

        i += 2 + natoms

    return np.array(frames), np.array(energies)


def _read_numfreq_result(numfreq_dir: Path) -> Dict[str, Any]:
    """Read NumFreq or AnFreq result, extract coords, energy, hessian, frequencies, thermochemistry

    Strategy:
    - Coords & Energy: Read from the OPTIMIZATION output (not from frequency calculation)
    - Hessian & Frequencies & Thermochemistry: Read from frequency calculation (NumFreq or AnFreq)
    - Gradient: Removed (not needed, optimization converged so gradient ≈ 0)

    Directory mapping:
    - 05_tsfreq (TS freq) → optimization in 04_tsopt → TS_opt.xyz
    - 07a_freq_irc_forward → optimization in 07_opt_irc_forward → IRC_forward_opt.xyz
    - 08a_freq_irc_backward → optimization in 08_opt_irc_backward → IRC_backward_opt.xyz
    """
    # Try NumFreq first, then AnFreq
    numfreq_result_file = numfreq_dir / "ASH_NumFreq.result"
    anfreq_result_file = numfreq_dir / "ASH_AnFreq.result"

    is_anfreq = False
    if numfreq_result_file.exists():
        result_file = numfreq_result_file
    elif anfreq_result_file.exists():
        result_file = anfreq_result_file
        is_anfreq = True
    else:
        return {}

    with open(result_file) as f:
        data = json.load(f)

    result = {
        "frequencies_cm-1": np.array(data.get("frequencies", [])),
        "thermochemistry": {
            "ZPVE_Eh": data["thermochemistry"]["ZPVE"],
            "Hcorr_Eh": data["thermochemistry"]["Hcorr"],
            "Gcorr_Eh": data["thermochemistry"]["Gcorr"],
        },
    }

    # Determine which optimization directory to read coords/energy from
    parent_dir = numfreq_dir.parent
    freq_dirname = numfreq_dir.name

    if freq_dirname == "05_tsfreq":
        opt_dir = parent_dir / "04_tsopt"
        opt_xyz = parent_dir / "TS_opt.xyz"
    elif freq_dirname == "07a_freq_irc_forward":
        opt_dir = parent_dir / "07_opt_irc_forward"
        opt_xyz = parent_dir / "IRC_forward_opt.xyz"
    elif freq_dirname == "08a_freq_irc_backward":
        opt_dir = parent_dir / "08_opt_irc_backward"
        opt_xyz = parent_dir / "IRC_backward_opt.xyz"
    else:
        # Fallback: try to infer from available files
        opt_dir = None
        opt_xyz = parent_dir / "TS_opt.xyz"
        if not opt_xyz.exists():
            opt_xyz = parent_dir / "IRC_forward_opt.xyz"
        if not opt_xyz.exists():
            opt_xyz = parent_dir / "IRC_backward_opt.xyz"

    # Read coordinates from optimized structure
    if opt_xyz.exists():
        _, coords = _read_xyz_elements_coords(opt_xyz)
        result["coords_A"] = coords

    # Read energy from optimization output
    if opt_dir and opt_dir.exists():
        # Try xTB output first
        xtb_out = opt_dir / "xtb_.out"
        if xtb_out.exists():
            with open(xtb_out) as f:
                for line in f:
                    if "TOTAL ENERGY" in line:
                        parts = line.split()
                        result["energy_Eh"] = float(parts[3])

        # Try ORCA output
        orca_out = opt_dir / "orca_.out"
        if orca_out.exists():
            with open(orca_out) as f:
                for line in f:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        result["energy_Eh"] = float(line.split()[-1])

    # Read Hessian from frequency calculation
    if is_anfreq:
        # AnFreq: Hessian is directly in numfreq_dir
        hess_file = numfreq_dir / "Hessian"
        if hess_file.exists():
            result["hessian_Eh_per_A2"] = np.loadtxt(hess_file)
    else:
        # NumFreq: Hessian is in Numfreq_dir subdirectory
        hess_file = numfreq_dir / "Numfreq_dir" / "Hessian"
        if not hess_file.exists():
            # Fallback: check directly in numfreq_dir
            hess_file = numfreq_dir / "Hessian"
        if hess_file.exists():
            result["hessian_Eh_per_A2"] = np.loadtxt(hess_file)

    return result


def _save_workflow_result_pkl(run_dir: Path, charge: int, mult: int, elements: List[str], temperature: float = 298.15, pressure: float = 1.0) -> Optional[Path]:
    """Collect all workflow results and save as pickle file.

    Returns None if any required data is missing.

    Note: All data (coords, energy, gradient, hessian, frequencies, thermochemistry)
    for each species comes from the same NumFreq calculation for consistency.
    """
    try:
        elements_r, r_coords = _read_xyz_elements_coords(run_dir / "R_opt.xyz")

        irc_traj, irc_energies = _read_irc_trajectory(run_dir / "06_irc/geometric_OPTtraj_irc.xyz")

        # Read NumFreq results (all data from the same frequency calculation)
        ts_numfreq = _read_numfreq_result(run_dir / "05_tsfreq")
        irc_f_numfreq = _read_numfreq_result(run_dir / "07a_freq_irc_forward")
        irc_b_numfreq = _read_numfreq_result(run_dir / "08a_freq_irc_backward")

        # Check for missing data (removed gradient_Eh_per_A from required keys)
        required_keys = ["coords_A", "energy_Eh", "hessian_Eh_per_A2",
                         "frequencies_cm-1", "thermochemistry"]
        for name, data in [("ts", ts_numfreq), ("irc_forward", irc_f_numfreq), ("irc_backward", irc_b_numfreq)]:
            if not data or any(k not in data for k in required_keys):
                print(f"[WARNING] Missing NumFreq data for {name}")
                return None

        if irc_traj is None or irc_energies is None:
            print(f"[WARNING] Missing IRC trajectory data")
            return None

        data = {
            "units": {
                "coords": "Å",
                "energy": "Eh",
                "hessian": "Eh/Å²",
                "frequencies": "cm⁻¹",
                "thermo": "Eh",
                "temperature": "K",
                "pressure": "atm",
            },
            "system": {
                "elements": elements_r,
                "charge": charge,
                "mult": mult,
                "natoms": len(elements_r),
                "temperature": temperature,
                "pressure": pressure,
            },
            # reactant = IRC backward endpoint (all from NumFreq)
            "reactant": {
                "coords_A": irc_b_numfreq["coords_A"],
                "energy_Eh": irc_b_numfreq["energy_Eh"],
                
                "hessian_Eh_per_A2": irc_b_numfreq["hessian_Eh_per_A2"],
                "frequencies_cm-1": irc_b_numfreq["frequencies_cm-1"],
                "thermochemistry": irc_b_numfreq["thermochemistry"],
            },
            # product = IRC forward endpoint (all from NumFreq)
            "product": {
                "coords_A": irc_f_numfreq["coords_A"],
                "energy_Eh": irc_f_numfreq["energy_Eh"],
                
                "hessian_Eh_per_A2": irc_f_numfreq["hessian_Eh_per_A2"],
                "frequencies_cm-1": irc_f_numfreq["frequencies_cm-1"],
                "thermochemistry": irc_f_numfreq["thermochemistry"],
            },
            # ts (all from NumFreq)
            "ts": {
                "coords_A": ts_numfreq["coords_A"],
                "energy_Eh": ts_numfreq["energy_Eh"],
                
                "hessian_Eh_per_A2": ts_numfreq["hessian_Eh_per_A2"],
                "frequencies_cm-1": ts_numfreq["frequencies_cm-1"],
                "thermochemistry": ts_numfreq["thermochemistry"],
            },
            "irc": {
                "coords_A": irc_traj,
                "energy_Eh": irc_energies,
            },
            # IRC forward endpoint (same as product, kept for backward compatibility)
            "irc_forward": {
                "coords_A": irc_f_numfreq["coords_A"],
                "energy_Eh": irc_f_numfreq["energy_Eh"],
                
                "hessian_Eh_per_A2": irc_f_numfreq["hessian_Eh_per_A2"],
                "frequencies_cm-1": irc_f_numfreq["frequencies_cm-1"],
                "thermochemistry": irc_f_numfreq["thermochemistry"],
            },
            # IRC backward endpoint (same as reactant, kept for backward compatibility)
            "irc_backward": {
                "coords_A": irc_b_numfreq["coords_A"],
                "energy_Eh": irc_b_numfreq["energy_Eh"],
                
                "hessian_Eh_per_A2": irc_b_numfreq["hessian_Eh_per_A2"],
                "frequencies_cm-1": irc_b_numfreq["frequencies_cm-1"],
                "thermochemistry": irc_b_numfreq["thermochemistry"],
            },
        }

        pkl_path = run_dir / "ash_xtb_result.pkl"
        with open(pkl_path, "wb") as f:
            pickle.dump(data, f)

        return pkl_path
    except Exception as e:
        print(f"[WARNING] Failed to save pkl result: {e}")
        return None


def _load_config(config_path: Path) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path) as f:
        return yaml.safe_load(f)


def _create_theory_from_config(theory_config: Dict[str, Any], theory_name: str = "theory") -> Any:
    """
    Factory function to create theory object based on config.

    Uses the modular theory_factories module for pluggable theory creation.

    Args:
        theory_config: Theory configuration dictionary
        theory_name: Name for logging (e.g., "main_theory", "neb_theory")

    Returns:
        Theory object (xTBTheory, ORCATheory, MACETheory, etc.)
    """
    from theory_factories import create_theory_from_config, list_theories

    theory = create_theory_from_config(theory_config)
    if theory is None:
        available = ", ".join(list_theories())
        raise ValueError(
            f"Failed to create theory with type '{theory_config.get('type')}'. "
            f"Available theories: {available}"
        )

    # Print info for backward compatibility
    theory_type = theory_config["type"]
    print(f"[INFO] Created {theory_name} with {theory_type}")
    print(f"       theorynamelabel: {theory.theorynamelabel}")

    return theory


def _get_freq_function_and_params(theory_type: str, freq_method: str, freq_npoint: int, freq_runmode: str, freq_cores: int, temperature: float, pressure: float):
    """
    Get appropriate frequency function and parameters based on theory type and freq_method.

    Args:
        theory_type: "xtb" or "orca"
        freq_method: "auto", "numerical", or "analytical"
        freq_npoint: npoint for numerical frequencies
        freq_runmode: runmode for numerical frequencies
        freq_cores: cores for numerical frequencies
        temperature: Temperature in K for thermochemistry
        pressure: Pressure in atm for thermochemistry

    Returns:
        Tuple of (freq_function, freq_kwargs)
    """
    # Determine which frequency method to use
    if freq_method == "auto":
        use_analytical = (theory_type == "orca")
    elif freq_method == "analytical":
        if theory_type != "orca":
            raise ValueError(f"Analytical frequencies only supported for ORCA, not {theory_type}")
        use_analytical = True
    elif freq_method == "numerical":
        use_analytical = False
    else:
        raise ValueError(f"Invalid freq_method: {freq_method}. Supported: auto, numerical, analytical")

    if use_analytical:
        from ash import AnFreq
        print(f"[INFO] Using AnFreq (analytical frequencies) for {theory_type}")
        print(f"       temp={temperature}K, pressure={pressure}atm")
        # AnFreq does NOT have numcores parameter - parallelization is handled by the theory object
        return AnFreq, {"temp": temperature, "pressure": pressure}
    else:
        from ash import NumFreq
        print(f"[INFO] Using NumFreq (numerical frequencies) for {theory_type}")
        print(f"       temp={temperature}K, pressure={pressure}atm")
        return NumFreq, {
            "npoint": freq_npoint,
            "runmode": freq_runmode,
            "numcores": freq_cores,
            "temp": temperature,
            "pressure": pressure,
        }


def _load_config_to_namespace(config: Dict[str, Any]) -> argparse.Namespace:
    """Convert config dictionary to argparse Namespace."""
    args = argparse.Namespace()

    # Map config structure to flat argument names
    args.reactant = config["inputs"]["reactant"]
    args.product = config["inputs"]["product"]
    args.ts_guess = config["inputs"]["ts_guess"]

    args.charge = config["system"]["charge"]
    args.mult = config["system"]["mult"]
    args.temperature = config["system"].get("temperature", 298.15)  # Default 298.15 K
    args.pressure = config["system"].get("pressure", 1.0)  # Default 1.0 atm

    # Theory configurations (store as dict for factory function)
    args.main_theory_config = config["main_theory"]
    args.neb_theory_config = config["neb_theory"]

    args.neb_images = config["neb"]["images"]
    args.neb_maxiter = config["neb"]["maxiter"]
    args.neb_interpolation = config["neb"]["interpolation"]
    args.neb_CI = config["neb"].get("CI", True)  # Default to True for backward compatibility
    args.neb_runmode = config["neb"]["runmode"]
    args.neb_cores = config["neb"]["cores"]
    args.neb_use_ts_guess = config["neb"].get("use_ts_guess", True)  # Default to True

    args.opt_maxiter = config["optimization"]["maxiter"]
    args.convergence_setting = config["optimization"].get("convergence_setting", "ORCA")

    args.tsopt_maxiter = config["tsopt"]["maxiter"]
    args.tsopt_convergence_setting = config["tsopt"].get("convergence_setting", "ORCA")
    args.tsopt_hessian = config["tsopt"]["hessian"]

    args.freq_method = config["frequency"].get("method", "auto")
    args.freq_npoint = config["frequency"]["npoint"]
    args.freq_runmode = config["frequency"]["runmode"]
    args.freq_cores = config["frequency"]["cores"]
    args.imag_threshold = config["frequency"]["imag_threshold"]

    args.irc_maxiter = config["irc"]["maxiter"]

    # Endpoint match configuration
    args.endpoint_match_config = config.get("endpoint_match", {"method": "smiles_openbabel"})

    args.outdir = config["output"]["base_dir"]
    args.printlevel = config["output"]["printlevel"]

    return args


def main() -> int:
    parser = argparse.ArgumentParser(
        description="ASH xTB explicit-step workflow smoke test (Opt R/P -> CI-NEB(+TS guess) -> TSOpt -> NumFreq -> IRC). "
                    "All parameters must be specified in a YAML config file."
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file (required)")

    cli_args = parser.parse_args()

    # Load config file
    repo_root, local_ash_root = _add_local_ash_to_syspath()
    config_path = (repo_root / cli_args.config).resolve() if not Path(cli_args.config).is_absolute() else Path(cli_args.config)

    if not config_path.exists():
        print(f"[ERROR] Config file not found: {config_path}")
        print(f"[ERROR] Please provide a valid config file using --config parameter")
        return 1

    config = _load_config(config_path)
    args = _load_config_to_namespace(config)

    from ash import Fragment, NEB, Optimizer  # noqa: E402
    # NumFreq and AnFreq will be imported dynamically based on theory type

    # Create two theory objects
    print("\n" + "="*80)
    print("Creating Theory Objects")
    print("="*80)

    main_theory = _create_theory_from_config(args.main_theory_config, "main_theory")
    neb_theory = _create_theory_from_config(args.neb_theory_config, "neb_theory")

    print("="*80 + "\n")

    in_reactant = (repo_root / args.reactant).resolve() if not Path(args.reactant).is_absolute() else Path(args.reactant)
    in_product = (repo_root / args.product).resolve() if not Path(args.product).is_absolute() else Path(args.product)
    in_ts_guess = (repo_root / args.ts_guess).resolve() if not Path(args.ts_guess).is_absolute() else Path(args.ts_guess)
    # Validate input files
    for p, label in [(in_reactant, "reactant"), (in_product, "product")]:
        if not p.exists():
            raise FileNotFoundError(f"Missing {label} xyz: {p}")

    # TS guess is optional - only validate if use_ts_guess=true
    if args.neb_use_ts_guess and not in_ts_guess.exists():
        raise FileNotFoundError(f"Missing ts_guess xyz: {in_ts_guess} (required when neb.use_ts_guess=true)")

    ts_tag = time.strftime("%Y%m%d_%H%M%S")
    run_dir = (repo_root / args.outdir / f"run_{ts_tag}").resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    summary: Dict[str, Any] = {
        "run_dir": str(run_dir),
        "timestamp": ts_tag,
        "inputs": {
            "reactant_xyz": str(in_reactant),
            "product_xyz": str(in_product),
            "ts_guess_xyz": str(in_ts_guess),
        },
        "params": vars(args),
        "local_ash_root": str(local_ash_root),
        "success": None,
        "stop_stage": None,
        "stop_reason": None,
        "exit_code": None,
        "stages": [],
    }

    stages: List[StageSummary] = []

    def _finalize(exit_code: int, stop_stage: Optional[str] = None, stop_reason: Optional[str] = None) -> int:
        summary["success"] = exit_code == 0
        summary["stop_stage"] = stop_stage
        summary["stop_reason"] = stop_reason
        summary["exit_code"] = exit_code
        out_path = _write_workflow_summary(run_dir, summary, stages)
        if exit_code == 0:
            print(f"[OK] Wrote workflow summary: {out_path}")
        else:
            print(f"[STOP] exit_code={exit_code} stage={stop_stage} reason={stop_reason}")
            print(f"[STOP] Wrote workflow summary: {out_path}")
        return exit_code

    # --- Stage 1: Opt reactant ---
    try:
        stage_dir = run_dir / "01_opt_reactant"
        with _pushd(stage_dir):
            frag_r = Fragment(xyzfile=str(in_reactant), charge=args.charge, mult=args.mult)
            res = Optimizer(
                fragment=frag_r,
                theory=main_theory,
                charge=args.charge,
                mult=args.mult,
                maxiter=args.opt_maxiter,
                printlevel=args.printlevel,
                convergence_setting=args.convergence_setting,
            )
            opt_xyz = stage_dir / "Fragment-optimized.xyz"
            _copy_if_exists(opt_xyz, run_dir / "R_opt.xyz")
            converged = _geometric_converged(stage_dir / "geometric_OPTtraj.log")
            stage_ok = converged and (run_dir / "R_opt.xyz").exists()
            file_count, files_txt = _run_stage_inventory(stage_dir, "opt_reactant")
            stages.append(
                StageSummary(
                    name="opt_reactant",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "energy": _as_float_or_none(getattr(res, "energy", None)),
                        "converged": converged,
                        "optimized_xyz": str(opt_xyz) if opt_xyz.exists() else None,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else "NOT_CONVERGED: reactant optimization did not converge",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="01_opt_reactant", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = _run_stage_inventory(stage_dir, "opt_reactant") if "stage_dir" in locals() else (0, run_dir / "missing.txt")
        stages.append(
            StageSummary(
                name="opt_reactant",
                workdir=str(stage_dir) if "stage_dir" in locals() else str(run_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="01_opt_reactant", stop_reason=stages[-1].error)

    # --- Stage 2: Opt product ---
    try:
        stage_dir = run_dir / "02_opt_product"
        with _pushd(stage_dir):
            frag_p = Fragment(xyzfile=str(in_product), charge=args.charge, mult=args.mult)
            res = Optimizer(
                fragment=frag_p,
                theory=main_theory,
                charge=args.charge,
                mult=args.mult,
                maxiter=args.opt_maxiter,
                printlevel=args.printlevel,
                convergence_setting=args.convergence_setting,
            )
            opt_xyz = stage_dir / "Fragment-optimized.xyz"
            _copy_if_exists(opt_xyz, run_dir / "P_opt.xyz")
            converged = _geometric_converged(stage_dir / "geometric_OPTtraj.log")
            stage_ok = converged and (run_dir / "P_opt.xyz").exists()
            file_count, files_txt = _run_stage_inventory(stage_dir, "opt_product")
            stages.append(
                StageSummary(
                    name="opt_product",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "energy": _as_float_or_none(getattr(res, "energy", None)),
                        "converged": converged,
                        "optimized_xyz": str(opt_xyz) if opt_xyz.exists() else None,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else "NOT_CONVERGED: product optimization did not converge",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="02_opt_product", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = _run_stage_inventory(stage_dir, "opt_product")
        stages.append(
            StageSummary(
                name="opt_product",
                workdir=str(stage_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="02_opt_product", stop_reason=stages[-1].error)

    # --- Stage 3: CI-NEB with optional TS guess insertion ---
    try:
        stage_dir = run_dir / "03_neb"
        with _pushd(stage_dir):
            r_opt_xyz = run_dir / "R_opt.xyz"
            p_opt_xyz = run_dir / "P_opt.xyz"
            if not r_opt_xyz.exists() or not p_opt_xyz.exists():
                raise FileNotFoundError("Missing R_opt.xyz or P_opt.xyz (previous stage did not produce expected output)")

            frag_r = Fragment(xyzfile=str(r_opt_xyz), charge=args.charge, mult=args.mult)
            frag_p = Fragment(xyzfile=str(p_opt_xyz), charge=args.charge, mult=args.mult)

            # Prepare TS guess file based on configuration
            ts_guess_file_param = None
            if args.neb_use_ts_guess:
                if not in_ts_guess.exists():
                    raise FileNotFoundError(f"TS guess file not found: {in_ts_guess} (required when neb.use_ts_guess=true)")
                ts_guess_file_param = str(in_ts_guess)

            res = NEB(
                reactant=frag_r,
                product=frag_p,
                theory=neb_theory,
                images=args.neb_images,
                CI=args.neb_CI,
                maxiter=args.neb_maxiter,
                interpolation=args.neb_interpolation,
                TS_guess_file=ts_guess_file_param,
                runmode=args.neb_runmode,
                numcores=args.neb_cores,
                charge=args.charge,
                mult=args.mult,
                printlevel=args.printlevel,
            )

            saddle = getattr(res, "saddlepoint_fragment", None)
            neb_result_json = _read_json_if_exists(stage_dir / "ASH_NEB.result") or {}
            neb_label = neb_result_json.get("label") or getattr(res, "label", None)
            neb_converged = saddle is not None and (neb_label is None or "fail" not in str(neb_label).lower())
            ts_guess_path = run_dir / "TS_guess_from_NEB.xyz"
            if neb_converged and saddle is not None:
                saddle.write_xyzfile(xyzfilename=str(ts_guess_path))

            file_count, files_txt = _run_stage_inventory(stage_dir, "neb")
            stages.append(
                StageSummary(
                    name="neb",
                    workdir=str(stage_dir),
                    ok=neb_converged,
                    result={
                        "energy": _as_float_or_none(getattr(res, "energy", None)),
                        "has_saddlepoint_fragment": saddle is not None,
                        "label": neb_label,
                        "mep_energies_dict": _serialize_mep_energies(getattr(res, "MEP_energies_dict", None)),
                        "ts_guess_xyz": str(ts_guess_path) if ts_guess_path.exists() else None,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if neb_converged else f"NOT_CONVERGED: NEB did not converge (label={neb_label!r})",
                )
            )
            if not neb_converged:
                return _finalize(exit_code=2, stop_stage="03_neb", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = _run_stage_inventory(stage_dir, "neb")
        stages.append(
            StageSummary(
                name="neb",
                workdir=str(stage_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="03_neb", stop_reason=stages[-1].error)

    # --- Stage 4: TS optimization ---
    try:
        stage_dir = run_dir / "04_tsopt"
        with _pushd(stage_dir):
            ts_guess = run_dir / "TS_guess_from_NEB.xyz"
            if not ts_guess.exists():
                raise FileNotFoundError("Missing TS_guess_from_NEB.xyz (NEB did not provide saddlepoint_fragment?)")

            frag_ts = Fragment(xyzfile=str(ts_guess), charge=args.charge, mult=args.mult)
            hess_opt = (args.tsopt_hessian or "").strip() or None
            res = Optimizer(
                fragment=frag_ts,
                theory=main_theory,
                charge=args.charge,
                mult=args.mult,
                TSOpt=True,
                hessian=hess_opt,
                maxiter=args.tsopt_maxiter,
                printlevel=args.printlevel,
                convergence_setting=args.tsopt_convergence_setting,
            )

            opt_xyz = stage_dir / "Fragment-optimized.xyz"
            _copy_if_exists(opt_xyz, run_dir / "TS_opt.xyz")
            converged = _geometric_converged(stage_dir / "geometric_OPTtraj.log")
            stage_ok = converged and (run_dir / "TS_opt.xyz").exists()
            file_count, files_txt = _run_stage_inventory(stage_dir, "tsopt")
            stages.append(
                StageSummary(
                    name="tsopt",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "energy": _as_float_or_none(getattr(res, "energy", None)),
                        "converged": converged,
                        "hessian_option": hess_opt,
                        "optimized_xyz": str(opt_xyz) if opt_xyz.exists() else None,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else "NOT_CONVERGED: TS optimization did not converge",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="04_tsopt", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = _run_stage_inventory(stage_dir, "tsopt")
        stages.append(
            StageSummary(
                name="tsopt",
                workdir=str(stage_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="04_tsopt", stop_reason=stages[-1].error)

    # --- Stage 5: NumFreq on TS ---
    freq_result = None
    try:
        stage_dir = run_dir / "05_tsfreq"
        with _pushd(stage_dir):
            ts_opt = run_dir / "TS_opt.xyz"
            if not ts_opt.exists():
                raise FileNotFoundError("Missing TS_opt.xyz (TSOpt did not produce expected output)")
            frag_ts = Fragment(xyzfile=str(ts_opt), charge=args.charge, mult=args.mult)

            # Get appropriate frequency function based on theory type
            FreqFunc, freq_kwargs = _get_freq_function_and_params(
                args.main_theory_config["type"],
                args.freq_method,
                args.freq_npoint,
                args.freq_runmode,
                args.freq_cores,
                args.temperature,
                args.pressure
            )

            freq_result = FreqFunc(
                fragment=frag_ts,
                theory=main_theory,
                printlevel=args.printlevel,
                **freq_kwargs
            )
            freqs_obj = getattr(freq_result, "frequencies", None)
            freqs = [float(f) for f in freqs_obj] if freqs_obj is not None else []
            imag_count, imag_freqs = _count_significant_imaginary(freqs, args.imag_threshold)

            # Copy Hessian file (handle both NumFreq and AnFreq directory structures)
            # NumFreq: stage_dir/Numfreq_dir/Hessian
            # AnFreq: stage_dir/Hessian (directly in stage_dir)
            hess_src = stage_dir / "Numfreq_dir" / "Hessian"
            if not hess_src.exists():
                hess_src = stage_dir / "Anfreq_dir" / "Hessian"
            if not hess_src.exists():
                hess_src = stage_dir / "Hessian"
            _copy_if_exists(hess_src, run_dir / "Hessian_NumFreq")
            wrote_hessian = (run_dir / "Hessian_NumFreq").exists()
            stage_ok = imag_count == 1 and wrote_hessian
            if not stage_ok:
                if imag_count != 1:
                    stage_error = (
                        f"FILTER_FAILED: TS must have exactly 1 imaginary freq < -{args.imag_threshold} cm^-1 (got {imag_count})"
                    )
                else:
                    stage_error = "OUTPUT_MISSING: Hessian_NumFreq not found (required for IRC)"
            else:
                stage_error = None
            file_count, files_txt = _run_stage_inventory(stage_dir, "numfreq")
            stages.append(
                StageSummary(
                    name="numfreq",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "imag_threshold_cm-1": float(args.imag_threshold),
                        "imaginary_freq_count_significant": imag_count,
                        "imaginary_frequencies_cm-1": imag_freqs[:10],
                        "first_imag_freq_cm-1": imag_freqs[0] if len(imag_freqs) == 1 else None,
                        "frequencies_cm-1_first20": freqs[:20],
                        "wrote_hessian_file": str(run_dir / "Hessian_NumFreq") if (run_dir / "Hessian_NumFreq").exists() else None,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=stage_error,
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="05_tsfreq", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = _run_stage_inventory(stage_dir, "numfreq")
        stages.append(
            StageSummary(
                name="numfreq",
                workdir=str(stage_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="05_tsfreq", stop_reason=stages[-1].error)

    # --- Stage 6: IRC (geomeTRIC) ---
    try:
        stage_dir = run_dir / "06_irc"
        with _pushd(stage_dir):
            ts_opt = run_dir / "TS_opt.xyz"
            if not ts_opt.exists():
                raise FileNotFoundError("Missing TS_opt.xyz (TSOpt did not produce expected output)")
            if freq_result is None or getattr(freq_result, "hessian", None) is None:
                raise RuntimeError("Missing frequency Hessian (required for IRC)")

            frag_ts = Fragment(xyzfile=str(ts_opt), charge=args.charge, mult=args.mult)
            res = Optimizer(
                fragment=frag_ts,
                theory=main_theory,
                charge=args.charge,
                mult=args.mult,
                irc=True,
                hessian=freq_result.hessian,
                maxiter=args.irc_maxiter,
                printlevel=args.printlevel,
            )
            irc_log = stage_dir / "geometric_OPTtraj.log"
            converged = _geometric_converged(irc_log)
            irc_log_meta = {
                "has_forward_converged_line": _file_contains_text(irc_log, "IRC forward direction converged"),
                "has_backward_start_line": _file_contains_text(irc_log, "IRC backward direction starts here"),
            }
            stage_ok = converged and irc_log_meta["has_forward_converged_line"] and irc_log_meta["has_backward_start_line"]

            endpoints_info = None
            if stage_ok:
                endpoints_info = _extract_geometric_irc_endpoints(
                    irc_log_path=irc_log,
                    irc_traj_xyz=stage_dir / "geometric_OPTtraj_irc.xyz",
                    out_forward_xyz=run_dir / "IRC_forward_end.xyz",
                    out_backward_xyz=run_dir / "IRC_backward_end.xyz",
                )

            file_count, files_txt = _run_stage_inventory(stage_dir, "irc")
            stages.append(
                StageSummary(
                    name="irc",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "energy": _as_float_or_none(getattr(res, "energy", None)),
                        "converged": converged,
                        "log_meta": irc_log_meta,
                        "endpoints": endpoints_info,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else "NOT_CONVERGED: IRC did not converge in both directions",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="06_irc", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = _run_stage_inventory(stage_dir, "irc")
        stages.append(
            StageSummary(
                name="irc",
                workdir=str(stage_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="06_irc", stop_reason=stages[-1].error)

    # --- Stage 7: Re-optimize IRC forward endpoint ---
    try:
        stage_dir = run_dir / "07_opt_irc_forward"
        with _pushd(stage_dir):
            end_xyz = run_dir / "IRC_forward_end.xyz"
            if not end_xyz.exists():
                raise FileNotFoundError("Missing IRC_forward_end.xyz (endpoint extraction failed)")
            frag = Fragment(xyzfile=str(end_xyz), charge=args.charge, mult=args.mult)
            res = Optimizer(
                fragment=frag,
                theory=main_theory,
                charge=args.charge,
                mult=args.mult,
                maxiter=args.opt_maxiter,
                printlevel=args.printlevel,
                convergence_setting=args.convergence_setting,
            )
            opt_xyz = stage_dir / "Fragment-optimized.xyz"
            _copy_if_exists(opt_xyz, run_dir / "IRC_forward_opt.xyz")
            converged = _geometric_converged(stage_dir / "geometric_OPTtraj.log")
            stage_ok = converged and (run_dir / "IRC_forward_opt.xyz").exists()
            file_count, files_txt = _run_stage_inventory(stage_dir, "opt_irc_forward")
            stages.append(
                StageSummary(
                    name="opt_irc_forward",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "energy": _as_float_or_none(getattr(res, "energy", None)),
                        "converged": converged,
                        "optimized_xyz": str(opt_xyz) if opt_xyz.exists() else None,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else "NOT_CONVERGED: IRC forward endpoint re-optimization did not converge",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="07_opt_irc_forward", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = (
            _run_stage_inventory(stage_dir, "opt_irc_forward") if "stage_dir" in locals() else (0, run_dir / "missing.txt")
        )
        stages.append(
            StageSummary(
                name="opt_irc_forward",
                workdir=str(stage_dir) if "stage_dir" in locals() else str(run_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="07_opt_irc_forward", stop_reason=stages[-1].error)

    # --- Stage 07a: Freq on optimized IRC forward endpoint (should have no imaginary freqs) ---
    try:
        stage_dir = run_dir / "07a_freq_irc_forward"
        with _pushd(stage_dir):
            end_opt_xyz = run_dir / "IRC_forward_opt.xyz"
            if not end_opt_xyz.exists():
                raise FileNotFoundError("Missing IRC_forward_opt.xyz (endpoint re-optimization did not produce expected output)")
            frag = Fragment(xyzfile=str(end_opt_xyz), charge=args.charge, mult=args.mult)

            # Get appropriate frequency function based on theory type
            FreqFunc, freq_kwargs = _get_freq_function_and_params(
                args.main_theory_config["type"],
                args.freq_method,
                args.freq_npoint,
                args.freq_runmode,
                args.freq_cores,
                args.temperature,
                args.pressure
            )

            freq_res = FreqFunc(
                fragment=frag,
                theory=main_theory,
                printlevel=args.printlevel,
                **freq_kwargs
            )
            freqs_obj = getattr(freq_res, "frequencies", None)
            freqs = [float(f) for f in freqs_obj] if freqs_obj is not None else []
            imag_any = [f for f in freqs if f < 0.0]
            imag_count_any = len(imag_any)
            stage_ok = imag_count_any == 0

            file_count, files_txt = _run_stage_inventory(stage_dir, "freq_irc_forward")
            stages.append(
                StageSummary(
                    name="freq_irc_forward",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "imaginary_freq_count": imag_count_any,
                        "imaginary_frequencies_cm-1": imag_any[:10],
                        "frequencies_cm-1_first20": freqs[:20],
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else f"FILTER_FAILED: IRC forward endpoint has {imag_count_any} imaginary frequencies",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="07a_freq_irc_forward", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = (
            _run_stage_inventory(stage_dir, "freq_irc_forward") if "stage_dir" in locals() else (0, run_dir / "missing.txt")
        )
        stages.append(
            StageSummary(
                name="freq_irc_forward",
                workdir=str(stage_dir) if "stage_dir" in locals() else str(run_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="07a_freq_irc_forward", stop_reason=stages[-1].error)

    # --- Stage 8: Re-optimize IRC backward endpoint ---
    try:
        stage_dir = run_dir / "08_opt_irc_backward"
        with _pushd(stage_dir):
            end_xyz = run_dir / "IRC_backward_end.xyz"
            if not end_xyz.exists():
                raise FileNotFoundError("Missing IRC_backward_end.xyz (endpoint extraction failed)")
            frag = Fragment(xyzfile=str(end_xyz), charge=args.charge, mult=args.mult)
            res = Optimizer(
                fragment=frag,
                theory=main_theory,
                charge=args.charge,
                mult=args.mult,
                maxiter=args.opt_maxiter,
                printlevel=args.printlevel,
                convergence_setting=args.convergence_setting,
            )
            opt_xyz = stage_dir / "Fragment-optimized.xyz"
            _copy_if_exists(opt_xyz, run_dir / "IRC_backward_opt.xyz")
            converged = _geometric_converged(stage_dir / "geometric_OPTtraj.log")
            stage_ok = converged and (run_dir / "IRC_backward_opt.xyz").exists()
            file_count, files_txt = _run_stage_inventory(stage_dir, "opt_irc_backward")
            stages.append(
                StageSummary(
                    name="opt_irc_backward",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "energy": _as_float_or_none(getattr(res, "energy", None)),
                        "converged": converged,
                        "optimized_xyz": str(opt_xyz) if opt_xyz.exists() else None,
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else "NOT_CONVERGED: IRC backward endpoint re-optimization did not converge",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="08_opt_irc_backward", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = (
            _run_stage_inventory(stage_dir, "opt_irc_backward") if "stage_dir" in locals() else (0, run_dir / "missing.txt")
        )
        stages.append(
            StageSummary(
                name="opt_irc_backward",
                workdir=str(stage_dir) if "stage_dir" in locals() else str(run_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="08_opt_irc_backward", stop_reason=stages[-1].error)

    # --- Stage 08a: Freq on optimized IRC backward endpoint (should have no imaginary freqs) ---
    try:
        stage_dir = run_dir / "08a_freq_irc_backward"
        with _pushd(stage_dir):
            end_opt_xyz = run_dir / "IRC_backward_opt.xyz"
            if not end_opt_xyz.exists():
                raise FileNotFoundError("Missing IRC_backward_opt.xyz (endpoint re-optimization did not produce expected output)")
            frag = Fragment(xyzfile=str(end_opt_xyz), charge=args.charge, mult=args.mult)

            # Get appropriate frequency function based on theory type
            FreqFunc, freq_kwargs = _get_freq_function_and_params(
                args.main_theory_config["type"],
                args.freq_method,
                args.freq_npoint,
                args.freq_runmode,
                args.freq_cores,
                args.temperature,
                args.pressure
            )

            freq_res = FreqFunc(
                fragment=frag,
                theory=main_theory,
                printlevel=args.printlevel,
                **freq_kwargs
            )
            freqs_obj = getattr(freq_res, "frequencies", None)
            freqs = [float(f) for f in freqs_obj] if freqs_obj is not None else []
            imag_any = [f for f in freqs if f < 0.0]
            imag_count_any = len(imag_any)
            stage_ok = imag_count_any == 0

            file_count, files_txt = _run_stage_inventory(stage_dir, "freq_irc_backward")
            stages.append(
                StageSummary(
                    name="freq_irc_backward",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result={
                        "imaginary_freq_count": imag_count_any,
                        "imaginary_frequencies_cm-1": imag_any[:10],
                        "frequencies_cm-1_first20": freqs[:20],
                    },
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None if stage_ok else f"FILTER_FAILED: IRC backward endpoint has {imag_count_any} imaginary frequencies",
                )
            )
            if not stage_ok:
                return _finalize(exit_code=2, stop_stage="08a_freq_irc_backward", stop_reason=stages[-1].error)
    except Exception as e:
        file_count, files_txt = (
            _run_stage_inventory(stage_dir, "freq_irc_backward") if "stage_dir" in locals() else (0, run_dir / "missing.txt")
        )
        stages.append(
            StageSummary(
                name="freq_irc_backward",
                workdir=str(stage_dir) if "stage_dir" in locals() else str(run_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )
        return _finalize(exit_code=2, stop_stage="08a_freq_irc_backward", stop_reason=stages[-1].error)

    # --- Stage 9: Endpoint match check (modular molecular matcher) ---
    try:
        stage_dir = run_dir / "09_endpoint_match"
        with _pushd(stage_dir):
            r_opt = run_dir / "R_opt.xyz"
            p_opt = run_dir / "P_opt.xyz"
            irc_f_opt = run_dir / "IRC_forward_opt.xyz"
            irc_b_opt = run_dir / "IRC_backward_opt.xyz"
            for pth in [r_opt, p_opt, irc_f_opt, irc_b_opt]:
                if not pth.exists():
                    raise FileNotFoundError(f"Missing required XYZ for endpoint match: {pth}")

            # Create matcher from config
            endpoint_matcher = _get_endpoint_matcher(args.endpoint_match_config)
            if endpoint_matcher is None:
                print(f"[WARNING] Unknown endpoint_match method: {args.endpoint_match_config.get('method')}, falling back to OpenBabel SMILES")
                from molecular_matchers import OpenBabelSmilesMatcher
                endpoint_matcher = OpenBabelSmilesMatcher()

            print(f"[INFO] Endpoint match method: {endpoint_matcher.name}")

            # Run cross-isomorphism check
            match_res = _cross_isomorphism(
                true_r=r_opt,
                true_p=p_opt,
                irc_r=irc_b_opt,
                irc_p=irc_f_opt,
                matcher=endpoint_matcher,
            )
            _write_text(stage_dir / "endpoint_match.json", json.dumps(match_res, indent=2, ensure_ascii=False))

            stage_ok = bool(match_res.get("success")) and match_res.get("endpoint_match") == "2-end match"
            file_count, files_txt = _run_stage_inventory(stage_dir, "endpoint_match")
            stages.append(
                StageSummary(
                    name="endpoint_match",
                    workdir=str(stage_dir),
                    ok=stage_ok,
                    result=match_res,
                    file_count=file_count,
                    files_txt=str(files_txt),
                    error=None
                    if stage_ok
                    else f"FILTER_FAILED: endpoint_match={match_res.get('endpoint_match')} (success={match_res.get('success')})",
                )
            )
    except Exception as e:
        file_count, files_txt = (
            _run_stage_inventory(stage_dir, "endpoint_match") if "stage_dir" in locals() else (0, run_dir / "missing.txt")
        )
        stages.append(
            StageSummary(
                name="endpoint_match",
                workdir=str(stage_dir) if "stage_dir" in locals() else str(run_dir),
                ok=False,
                result={},
                file_count=file_count,
                files_txt=str(files_txt),
                error=f"EXCEPTION: {e}",
            )
        )

    # --- Save workflow result as PKL ---
    elements_r, _ = _read_xyz_elements_coords(run_dir / "R_opt.xyz")
    pkl_path = _save_workflow_result_pkl(run_dir, charge=args.charge, mult=args.mult, elements=elements_r, temperature=args.temperature, pressure=args.pressure)
    if pkl_path:
        print(f"[OK] Saved workflow result PKL: {pkl_path}")
    else:
        print("[WARNING] Failed to save workflow result PKL (missing data)")

    return _finalize(exit_code=0)


if __name__ == "__main__":
    raise SystemExit(main())
