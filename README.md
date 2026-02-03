# Elementary-Reaction-Validator

Automated workflow for chemical reaction pathway validation using ASH (A Suite for High-throughput quantum chemistry).

## Features

- **Multiple Theory Levels**: xTB, ORCA, MACE, FairChem, DPA
- **Complete 9-Stage Workflow**: From reactant/product optimization to IRC endpoint validation
- **Modular Endpoint Matching**: SMILES, Graph Isomorphism, RMSD, SOAP
- **Configurable Thermochemistry**: Adjustable temperature and pressure

## Quick Start

```bash
# Activate environment (choose based on theory)
conda activate xxx

python ash_workflow.py --config config.yaml
```

## Environment Setup

### Core Installation (Required)
```bash
pip install git+https://github.com/RagnarB83/ash.git@NEW
pip install pyyaml

conda install -c conda-forge openbabel
pip install rmsd dscribe networkx pymatgen rdkit
```

### Theory-Specific

| Theory | Python | Installation |
|--------|--------|--------------|
| **xTB** | 3.10 | `pip install xtb-python` |
| **ORCA** | 3.10 | Install manually: http://sobereva.com/451 |
| **MACE** | 3.9 | `pip install mace-torch` |
| **UMA (FairChem)** | 3.10 | `pip install fairchem-core` |
| **DPA** | 3.12 | `pip install deepmd-kit[torch]` |

## 9-Stage Workflow

| Stage | Description |
|-------|-------------|
| 1. opt_reactant | Optimize reactant geometry |
| 2. opt_product | Optimize product geometry |
| 3. neb | NEB/CI-NEB path search, generate TS guess |
| 4. tsopt | Transition state optimization |
| 5. numfreq | TS frequency (confirm 1 imaginary mode) |
| 6. irc | IRC calculation from TS |
| 7. opt_irc_forward + freq | Optimize and validate forward endpoint |
| 8. opt_irc_backward + freq | Optimize and validate backward endpoint |
| 9. endpoint_match | Cross-validate endpoints vs reactant/product |

## Configuration

All parameters are specified in the YAML config file (`ash_xtb_workflow_config.yaml`).

## Output Results

Run directory: `ash_smoketest_runs/run_YYYYMMDD_HHMMSS/`

### Key Output Files

| File | Description |
|------|-------------|
| `workflow_summary.json` | Human-readable summary with all stage results |
| `ash_result.pkl` | Binary data with NumPy arrays |
| `TS_guess_from_NEB.xyz` | TS geometry from NEB |
| `IRC_forward_end.xyz` | Forward IRC endpoint |
| `IRC_backward_end.xyz` | Backward IRC endpoint |


## Version

Current version: **3.4** (2026-02-03)
