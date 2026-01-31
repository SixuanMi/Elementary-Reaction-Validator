# Elementary-Reaction-Validator

Automated workflow for chemical reaction pathway validation using ASH (A Suite for High-throughput quantum chemistry).

## Features

- **Multiple Theory Levels**: xTB, ORCA, MACE, FairChem, DPA - with independent settings for main and NEB calculations
- **Flexible Frequency Methods**: Automatic selection of NumFreq (numerical) or AnFreq (analytical) based on theory
- **Complete 9-Stage Workflow**: From reactant/product optimization to IRC endpoint validation
- **Modular Endpoint Matching**: 5 different molecular matching methods (SMILES, Graph Isomorphism, RMSD, SOAP)
- **Configurable Thermochemistry**: Adjustable temperature and pressure for free energy corrections
- **Comprehensive Output**: JSON summary + PKL binary data with all results

## Quick Start

```bash
# Activate the appropriate environment for your theory
conda activate dpa3        # For DPA
# conda activate mace      # For MACE/FairChem
# conda activate ash-py310 # For xTB/ORCA

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
  type: "xtb"  # xtb, orca, mace, fairchem, or dpa
  numcores: 1
  printlevel: 1

  # xTB settings (only used when type=xtb)
  xtb:
    method: "GFN2"        # GFN1, GFN2
    runmode: "inputfile"  # inputfile or library

  # ORCA settings (only used when type=orca)
  orca:
    orcasimpleinput: "! B97-3c"
    orcablocks: |
      %scf
        maxiter 200
      end
    moreadfile: null  # Optional: path to .gbw file for orbital guess

  # MACE settings (only used when type=mace)
  mace:
    model_file: "path/to/mace-mh-1.model"  # Path to MACE model file
    device: "cpu"                            # cpu, cuda, opencl, mps

  # FairChem settings (only used when type=fairchem)
  fairchem:
    model_file: "path/to/uma-s-1p1.pt"      # Path to FairChem model file
    task_name: "omol"                        # oc20, omol, omat, odac, omc
    device: "cpu"                            # cpu, cuda, opencl, mps

  # DPA3 settings (only used when type=dpa)
  dpa:
    model_file: "path/to/DPA-3.1-3M.pt"      # Path to DPA model file
    head: "Organic_Reactions"                # Model head for multitask models
    device: "cpu"                            # cpu, cuda, opencl, mps

# NEB theory (can be cheaper/faster than main_theory)
neb_theory:
  type: "xtb"
  numcores: 1
  # ... same theory-specific settings as main_theory ...

# NEB settings
neb:
  images: 8
  maxiter: 200
  interpolation: "IDPP"
  CI: true
  runmode: "serial"
  cores: 1

# Optimization settings
optimization:
  maxiter: 200
  convergence_setting: "ORCA"

# TS optimization settings
tsopt:
  maxiter: 200
  convergence_setting: "ORCA"
  hessian: null  # xtb: "xtb"; ORCA: "numerical" or analytical; MACE/FairChem: null

# Frequency settings
frequency:
  method: "auto"  # auto (ORCA→AnFreq, xTB/ML→NumFreq), numerical, analytical
  npoint: 2
  runmode: "serial"
  cores: 1
  imag_threshold: 20.0  # cm^-1

# IRC settings
irc:
  maxiter: 200

# Endpoint matching
endpoint_match:
  method: "smiles_openbabel"  # smiles_openbabel, smiles_rdkit, graph_isomorphism, rmsd, soap
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

## Theory Types

### xTB (GFN)
Semi-empirical tight-binding method. Fast but less accurate.

```yaml
main_theory:
  type: "xtb"
  xtb:
    method: "GFN2"        # GFN1 or GFN2
    runmode: "inputfile"  # inputfile (disk) or library (memory)
```

**Environment**: `ash-py310` or any with ASH installed

### ORCA
DFT with various functionals. Slower but more accurate.

```yaml
main_theory:
  type: "orca"
  orca:
    orcasimpleinput: "! B97-3c"     # Method and basis
    orcablocks: |
      %scf
        maxiter 200
      end
    moreadfile: "previous.gbw"     # Optional orbital guess
```

**Environment**: `ash-py310` with ORCA installed

### MACE
Machine learning potential for molecules. Fast and accurate for organic molecules.

```yaml
main_theory:
  type: "mace"
  mace:
    model_file: "mlps/mace-mh-1.model"  # Path to .model file
    device: "cpu"                        # cpu, cuda, opencl, mps
```

**Environment**: `mace`

**Installation**:
```bash
pip install mace-torch
# Optional: pip install cuequivariance_torch  # Acceleration
```

**Models**:
- MACE-OFF24: General-purpose organic molecules
- MACE-MP: Materials and periodic systems
- Custom models: Train your own with MACE

### FairChem (FAIR-Chem)
Universal ML potential for catalysis and organic chemistry.

```yaml
main_theory:
  type: "fairchem"
  fairchem:
    model_file: "mlps/uma-s-1p1.pt"  # Path to .pt file
    task_name: "omol"                # oc20, omol, omat, odac, omc
    device: "cpu"                    # cpu, cuda, opencl, mps
```

**Environment**: `mace` (shared dependencies with MACE)

**Task Names**:
- `oc20`: OC20 (catalysis)
- `omol`: OC22 molecular (organic molecules)
- `omat`: OC22 materials
- `odac`: DAC (catalysis)
- `omc`: MC (materials)

**Installation**:
```bash
pip install fair-chem
```

### DPA3 (Deep Potential Approximation)
Universal ML potential with multitask support for diverse chemical systems.

```yaml
main_theory:
  type: "dpa"
  dpa:
    model_file: "mlps/DPA-3.1-3M.pt"  # Path to .pt file
    head: "Organic_Reactions"         # Model head for multitask models
    device: "cpu"                     # cpu, cuda, opencl, mps
```

**Environment**: `dpa3`

**Available Heads** (for DPA-3.1-3M multitask model):
- `Organic_Reactions`: Organic molecular reactions
- `OC20M`: OC20 catalysis
- `OC22`: OC22 molecular and materials
- `Omat24`: Materials
- `SPICE2`: Molecular datasets
- And many more (use `dp --pt show model.pt model-branch` to list)

**Installation**:
```bash
# DPA3 is part of DeepMD-kit
conda install -c deepmodeling deepmd
```

### Theory Comparison

| Theory | Speed | Accuracy | Use Case | Environment |
|--------|-------|----------|----------|-------------|
| xTB GFN2 | Fast | Medium | Screening, large systems | ash-py310 |
| ORCA DFT | Slow | High | Benchmarking, final results | ash-py310 |
| MACE | Fast | High | Organic molecules, reactions | mace |
| FairChem | Medium | High | Catalysis, diverse chemistry | mace |
| DPA3 | Medium | High | Multitask, diverse systems | dpa3 |

### Mixed Theory Strategies

Use cheaper theory for NEB, accurate theory for final TS:

```yaml
main_theory:
  type: "orca"      # Accurate for opt, tsopt, freq, irc

neb_theory:
  type: "mace"      # Fast for path search
```

Or use ML for exploration, DFT for validation:

```yaml
main_theory:
  type: "mace"      # Fast screening

neb_theory:
  type: "xtb"       # Even faster initial path
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
  "timestamp": "20260131_065434",
  "success": true,
  "stages": [
    {
      "name": "opt_reactant",
      "ok": true,
      "result": {
        "energy": -1.529,
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
        "method": "smiles_openbabel",
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
  "method": "smiles_openbabel",
  "smi_true_r": "CC=O",
  "smi_true_p": "OC=C",
  "smi_irc_r": "CC=O",
  "smi_irc_p": "OC=C"
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

## Theory Factory Architecture

The modular theory system uses a Factory Pattern:

```
theory_factories/
├── __init__.py              # Main exports
├── base.py                  # TheoryFactoryBase abstract class
├── registry.py              # Theory registration and factory
└── theories/
    ├── semiempirical/
    │   └── xtb_factory.py
    ├── dft/
    │   └── orca_factory.py
    └── ml/
        ├── mace_factory.py
        ├── fairchem_factory.py
        └── dpa_factory.py
```

## Dependencies

### Core (Required)
- Python >= 3.9
- ASH (A Suite for High-throughput quantum chemistry)
- NumPy, PyYAML

### Theory-Specific

**xTB**: Included with ASH

**ORCA**: Must be installed separately

**MACE**:
```bash
pip install mace-torch
# Optional: pip install cuequivariance_torch
```

**FairChem**:
```bash
pip install fair-chem
```

**DPA3**:
```bash
conda install -c deepmodeling deepmd
```

### Molecular Matchers (Optional - install as needed)

All matcher dependencies can be installed with pip:

```bash
pip install rmsd dscribe networkx pymatgen rdkit
```

For OpenBabel SMILES matcher, use conda:

```bash
conda install -c conda-forge openbabel
```

### Recommended Environments

**For xTB/ORCA**:
```bash
conda create -n ash-py310 python=3.10
conda activate ash-py310
# Install ASH and matcher dependencies
pip install rmsd dscribe networkx pymatgen rdkit
conda install -c conda-forge openbabel
```

**For MACE/FairChem**:
```bash
conda create -n mace python=3.9
conda activate mace
# Install ASH
pip install mace-torch fair-chem
# Optional: pip install cuequivariance_torch
```

**For DPA3**:
```bash
conda create -n dpa3 python=3.12
conda activate dpa3
# Install ASH
conda install -c deepmodeling deepmd
```

## Version

Current version: **3.3** (2026-01-31)

### Recent Changes (v3.3)
- DPA3 theory support via modular theory factory
- Multitask model support with configurable heads
- Multi-environment support (dpa3 for DPA, mace for MACE/FairChem, ash-py310 for xTB/ORCA)

### Previous Versions
- v3.2: MACE theory support, FairChem theory support
- v3.1: Modular molecular matcher system with 5 methods
- v3.0: Dual theory levels, analytical frequency support
- v2.0: Configuration file only, parameter modularization
- v1.0: Basic xTB workflow
