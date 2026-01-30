#!/usr/bin/env python3
"""Test script to verify config loading works correctly."""

from pathlib import Path
import sys

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML is not installed.")
    print("Please install it with: pip install pyyaml")
    print("Or in conda environment: conda install pyyaml")
    sys.exit(1)

config_path = Path("ash_xtb_workflow_config.yaml")

if not config_path.exists():
    print(f"ERROR: Config file not found: {config_path}")
    sys.exit(1)

with open(config_path) as f:
    config = yaml.safe_load(f)

print("✓ Config file loaded successfully!")
print("\nConfiguration summary:")
print(f"  Reactant: {config['inputs']['reactant']}")
print(f"  Product: {config['inputs']['product']}")
print(f"  TS guess: {config['inputs']['ts_guess']}")
print(f"  xTB method: {config['xtb']['method']}")
print(f"  NEB images: {config['neb']['images']}")
print(f"  NEB maxiter: {config['neb']['maxiter']}")
print(f"  Optimization maxiter: {config['optimization']['maxiter']}")
print(f"  Frequency npoint: {config['frequency']['npoint']}")
print(f"  Output directory: {config['output']['base_dir']}")
print("\n✓ All config sections are valid!")
