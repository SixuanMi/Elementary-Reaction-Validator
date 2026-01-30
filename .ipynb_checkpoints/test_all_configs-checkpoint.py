#!/usr/bin/env python3
"""Quick test to verify the smoketest script can load config correctly."""

import sys
from pathlib import Path

# Add local ash to path (same as smoketest)
repo_root = Path(__file__).resolve().parent
local_ash_root = repo_root / "ash"
sys.path.insert(0, str(local_ash_root))

import yaml
import argparse

def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)

def test_config(config_name):
    print(f"\n{'='*60}")
    print(f"Testing: {config_name}")
    print('='*60)

    config_path = repo_root / config_name
    if not config_path.exists():
        print(f"❌ Config file not found: {config_path}")
        return False

    try:
        config = load_config(config_path)
        print("✓ Config loaded successfully")

        # Display key parameters
        print(f"\nKey parameters:")
        print(f"  NEB images: {config['neb']['images']}")
        print(f"  NEB maxiter: {config['neb']['maxiter']}")
        print(f"  Opt maxiter: {config['optimization']['maxiter']}")
        print(f"  TS maxiter: {config['tsopt']['maxiter']}")
        print(f"  Freq npoint: {config['frequency']['npoint']}")
        print(f"  IRC maxiter: {config['irc']['maxiter']}")
        print(f"  xTB cores: {config['xtb']['cores']}")

        return True
    except Exception as e:
        print(f"❌ Error loading config: {e}")
        return False

if __name__ == "__main__":
    configs = [
        "ash_xtb_workflow_config.yaml",
        "ash_xtb_workflow_config_fast.yaml",
        "ash_xtb_workflow_config_accurate.yaml",
    ]

    print("Testing all configuration files...")
    results = []

    for config in configs:
        results.append(test_config(config))

    print(f"\n{'='*60}")
    print("Summary:")
    print('='*60)
    for config, result in zip(configs, results):
        status = "✓ PASS" if result else "❌ FAIL"
        print(f"{status}: {config}")

    if all(results):
        print("\n✓ All configuration files are valid!")
        sys.exit(0)
    else:
        print("\n❌ Some configuration files have errors")
        sys.exit(1)
