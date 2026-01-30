# Elementary-Reaction-Validator

Automated workflow for chemical reaction pathway validation using ASH (A Suite for High-throughput quantum chemistry).

## Features

- **Dual Theory Levels**: Independent theory settings for main calculations and NEB path search
- **Flexible Frequency Methods**: Automatic selection of NumFreq (numerical) or AnFreq (analytical) based on theory
- **Complete 9-Stage Workflow**: From reactant/product optimization to IRC endpoint validation
- **Modular Endpoint Matching**: 5 different molecular matching methods (SMILES, Graph Isomorphism, RMSD, SOAP)
- **Configurable Thermochemistry**: Adjustable temperature and pressure for free energy corrections
- **Comprehensive Output**: JSON summary + PKL binary data with all results

## Quick Start

```bash
python ash_xtb_explicit_workflow_smoketest.py --config ash_xtb_workflow_config.yaml
```

## Configuration

All parameters are specified in the YAML config file (`ash_xtb_workflow_config.yaml`).

### Basic Structure

```yaml
# Input files
inputs:
  reactant: "xyz/reactant.xyz"
  product: "xyz/product.xyz"
  ts_guess: "xyz/ts_ini.xyz"

# System properties
system:
  charge: 0
  mult: 1
  temperature: 298.15
  pressure: 1.0

# Main theory (for opt, tsopt, freq, irc)
main_theory:
  type: "xtb"  # or "orca"
  numcores: 1
  xtb:
    method: "GFN2"
    runmode: "inputfile"

# NEB theory (can be cheaper/faster)
neb_theory:
  type: "xtb"
  numcores: 1
  xtb:
    method: "GFN2"
    runmode: "inputfile"

# NEB settings
neb:
  images: 8
  maxiter: 200
  interpolation: "IDPP"
  CI: true

# Endpoint matching
endpoint_match:
  method: "smiles_rdkit"  # Choose method
  smiles_rdkit:
    atom_map: false
    all_hs_explicit: true
  rmsd:
    threshold: 0.5
    allow_reorder: false
    heavy_only: false
  soap:
    r_cut: 6.0
    n_max: 8
    l_max: 6
    sigma: 1.0
    kernel_metric: "linear"
    threshold_similarity: 0.99

# Output
output:
  base_dir: "ash_xtb_smoketest_runs"
  printlevel: 1
```

## Endpoint Matching Methods

Stage 9 validates that IRC endpoints match the original reactant/product. Multiple methods are available:

| Method | Description | Strictness | Use Case |
|--------|-------------|------------|----------|
| `smiles_openbabel` | OpenBabel SMILES comparison | High | General use, stereochemistry |
| `smiles_rdkit` | RDKit SMILES with bond strategies | High | Complex molecules, atom mapping |
| `graph_isomorphism` | Pymatgen + NetworkX graph matching | High | Topology changes |
| `rmsd` | 3D RMSD with Kabsch alignment | Very High | Geometric comparison |
| `soap` | SOAP descriptor similarity | Medium | Similar conformations |

### Method Details

#### SMILES (OpenBabel)
```yaml
endpoint_match:
  method: "smiles_openbabel"
```
- Always includes stereochemistry
- Fast and reliable for most organic molecules

#### SMILES (RDKit)
```yaml
endpoint_match:
  method: "smiles_rdkit"
  smiles_rdkit:
    atom_map: false        # Include atom map numbers (AAM)
    all_hs_explicit: true  # true=explicit H, false=implicit H
```
- Multiple bond determination strategies (Hueckel, Vdw, Basic)
- Robust for complex molecules
- Optional atom mapping for reaction tracking

#### Graph Isomorphism
```yaml
endpoint_match:
  method: "graph_isomorphism"
  graph_isomorphism:
    local_env: "OpenBabelNN"  # or CovalentBondNN, MinimumDistanceNN
```
- Based on molecular graph topology
- Uses Pymatgen MoleculeGraph + NetworkX
- Good for detecting bond changes

#### RMSD
```yaml
endpoint_match:
  method: "rmsd"
  rmsd:
    threshold: 0.5         # Angstrom
    allow_reorder: false   # Hungarian algorithm for atom correspondence
    heavy_only: false      # Ignore hydrogen atoms
```
- 3D geometric comparison after Kabsch alignment
- Very strict - suitable for distinguishing similar structures
- Higher RMSD = more different structures

#### SOAP
```yaml
endpoint_match:
  method: "soap"
  soap:
    r_cut: 6.0                      # Cutoff radius (Angstrom)
    n_max: 8                        # Radial basis functions
    l_max: 6                        # Angular momentum
    sigma: 1.0                      # Gaussian width (higher = smoother)
    kernel_metric: "linear"         # "linear" or "rbf"
    kernel_gamma: null              # RBF gamma (null=auto)
    threshold_similarity: 0.99      # Similarity threshold
```
- Based on local atomic environments
- Uses dscribe SOAP + AverageKernel
- Lower sigma = more sensitive to small differences
- Suitable for conformational comparison

### Match Categories

| Category | Condition | Meaning |
|----------|-----------|---------|
| `2-end match` | true_r≟irc_r AND true_p≟irc_p | Correct IRC direction |
| `1-end match` | Partial match | Possible swapped/incorrect IRC |
| `No match` | No matches | IRC failed or wrong TS |

### Reaction Status

- `Chemical reaction`: IRC endpoints are chemically different
- `Conformational change`: IRC endpoints are similar (same connectivity)

## Output Files

### Main Output Directory
```
ash_xtb_smoketest_runs/run_YYYYMMDD_HHMMSS/
├── 01_opt_reactant/          # Reactant optimization
├── 02_opt_product/           # Product optimization
├── 03_neb/                   # NEB path search
├── 04_tsopt/                 # TS optimization
├── 05_tsfreq/                # TS frequency (imaginary mode check)
├── 06_irc/                   # IRC calculation
├── 07_opt_irc_forward/       # Forward endpoint optimization
├── 07a_freq_irc_forward/     # Forward endpoint frequency
├── 08_opt_irc_backward/      # Backward endpoint optimization
├── 08a_freq_irc_backward/    # Backward endpoint frequency
├── 09_endpoint_match/        # Endpoint validation
├── workflow_summary.json     # Human-readable summary
├── ash_xtb_result.pkl        # Binary data with NumPy arrays
└── [various XYZ files]
```

### workflow_summary.json

Contains complete results for each stage:

```json
{
  "run_dir": "...",
  "timestamp": "20260130_140506",
  "success": true,
  "stages": [
    {
      "name": "opt_reactant",
      "ok": true,
      "result": {
        "energy": -21.56673926153,
        "converged": true,
        "optimized_xyz": "..."
      }
    },
    ...
    {
      "name": "endpoint_match",
      "ok": true,
      "result": {
        "endpoint_match": "2-end match",
        "matches": [true, true, false, false],
        "rxn_status": "Chemical reaction",
        "method": "smiles_rdkit",
        "metadata": {...}
      }
    }
  ]
}
```

### endpoint_match.json

Located in `09_endpoint_match/`, contains detailed match results:

```json
{
  "performed": true,
  "success": true,
  "endpoint_match": "2-end match",
  "matches": [true, true, false, false],
  "rxn_status": "Chemical reaction",
  "method": "smiles_rdkit",
  "smi_true_r": "CC(=O)O",
  "smi_true_p": "CC(O)=O",
  "smi_irc_r": "CC(O)=O",
  "smi_irc_p": "CC(=O)O"
}
```

## Molecular Matcher Architecture

The modular matcher system uses a Strategy Pattern:

```
molecular_matchers/
├── __init__.py           # Main exports
├── base.py               # MoleculeMatcher abstract base class
├── registry.py           # Matcher registration and factory
├── smiles/
│   ├── openbabel_matcher.py
│   └── rdkit_matcher.py
├── graph/
│   └── isomorphism_matcher.py
├── geometry/
│   └── rmsd_matcher.py
└── descriptor/
    └── soap_matcher.py
```

### Creating a Custom Matcher

```python
from molecular_matchers import MoleculeMatcher, register_matcher

@register_matcher("my_method")
class MyMatcher(MoleculeMatcher):
    def compute_signature(self, xyz_file):
        # Compute signature from XYZ
        return signature

    def compare(self, sig1, sig2):
        # Compare two signatures
        return similarity >= threshold
```

## Dependencies

### Core (Required)
- Python >= 3.10
- ASH (A Suite for High-throughput quantum chemistry)
- NumPy, PyYAML

### Molecular Matchers (Optional - install as needed)

All matcher dependencies can be installed with pip:

```bash
pip install rmsd dscribe networkx pymatgen rdkit
```

For OpenBabel SMILES matcher, use conda:

```bash
conda install -c conda-forge openbabel
```

### Recommended Environment
```bash
conda create -n ash-py310 python=3.10
conda activate ash-py310
# Install ASH (per your setup)
# Then install matcher dependencies:
pip install rmsd dscribe networkx pymatgen rdkit
conda install -c conda-forge openbabel
```

## Version

Current version: **3.1** (2026-01-30)

### Recent Changes (v3.1)
- Modular molecular matcher system with 5 methods
- SOAP matcher using dscribe AverageKernel
- Graph isomorphism using Pymatgen + NetworkX
- RDKit SMILES with atom mapping support
- RMSD with Kabsch alignment and Hungarian reordering

### Previous Versions
- v3.0: Dual theory levels, analytical frequency support
- v2.0: Configuration file only, parameter modularization
- v1.0: Basic xTB workflow
