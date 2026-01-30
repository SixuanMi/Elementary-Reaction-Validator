# Elementary-Reaction-Validator

Automated workflow for chemical reaction pathway validation using ASH (A Suite for High-throughput quantum chemistry).

## Features

- **Dual Theory Levels**: Independent theory settings for main calculations and NEB path search
- **Flexible Frequency Methods**: Automatic selection of NumFreq (numerical) or AnFreq (analytical) based on theory
- **Complete 9-Stage Workflow**: From reactant/product optimization to IRC endpoint validation
- **Configurable Thermochemistry**: Adjustable temperature and pressure for free energy corrections
- **Comprehensive Output**: JSON summary + PKL binary data with all results

## Quick Start

```bash
python ash_xtb_explicit_workflow_smoketest.py --config ash_xtb_workflow_config.yaml
```

## Configuration

All parameters are specified in the YAML config file. See `CONFIG_README.md` for details.

Key sections:
- `inputs`: XYZ files for reactant, product, TS guess
- `system`: Charge, multiplicity, temperature, pressure
- `main_theory`: ORCA or xTB settings for optimization/freq/IRC
- `neb_theory`: Theory for NEB path search (can be cheaper)
- `frequency`: Method selection (auto/numerical/analytical)

## Output Files

- `workflow_summary.json`: Human-readable summary with all energies and diagnostics
- `ash_xtb_result.pkl`: Binary file with NumPy arrays (coords, Hessian, frequencies, etc.)
- `Hessian_NumFreq`: Hessian matrix for IRC calculation
- Stage directories (`01_opt_reactant`, `02_opt_product`, etc.)

## Documentation

- `ASH_XTB_WORKFLOW_NOTES.md`: Detailed workflow description and technical notes
- `CONFIG_README.md`: Configuration file format and parameters
- `CHANGELOG.md`: Version history and changes

## Version

Current version: **3.0** (2026-01-30)

See `CHANGELOG.md` for recent changes.
