# HamDiv Package Installation Guide

## Installation with uv

You can install this package using `uv` in several ways:

### Install from local directory (development)

If you're in the HamDiv directory:
**Manual:**
```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

This installs the package in "editable" mode, meaning changes to the source code will be reflected immediately.

### Install dependencies only

To install just the dependencies without the package:
```bash
uv pip install -r requirements.txt
```

### Install with development dependencies

For development work:
```bash
uv pip install -e ".[dev]"
```

## Usage

After installation, you can use the package:

```python
from hamdiv import diversity_all, HamDiv, diversity_vector, diversity_vector_array

# Example SMILES
smiles = [
    'Cc1cc(C(O)CNC2(C)CC2S)ccc1O',
    'CCc1cc(C(O)CNC(C)(C)C)ccc1O',
    'CCCC(C)NCC(O)c1ccc(O)c(CO)c1',
    'CNCC(S)c1ccc(O)c(CO)c1',
    'CNCC(C)(C)C(O)c1ccc(S)c(CO)c1'
]

# Calculate individual diversity metrics
richness = diversity_all(smiles=smiles, mode="Richness")
ham_div = diversity_all(smiles=smiles, mode="HamDiv")
# or directly:
ham_div_direct = HamDiv(smiles=smiles)

# Get all metrics as a vector (NEW!)
diversity_dict = diversity_vector(smiles=smiles)
metric_names, values_array = diversity_vector_array(smiles=smiles)
```

## Available Diversity Metrics

- `Richness`: Number of unique molecules
- `IntDiv`: Internal diversity 
- `SumDiv`: Sum diversity
- `Diam`: Diameter
- `SumDiam`: Sum diameter
- `Bot`: Bottleneck
- `SumBot`: Sum bottleneck
- `DPP`: Determinantal point process
- `NCircles-X`: N-circles with threshold X
- `HamDiv`: Hamiltonian diversity (main contribution)
- `FG`: Functional groups diversity
- `RS`: Ring systems diversity
- `BM`: Bemis-Murcko scaffolds diversity

## Dependencies

- numpy
- rdkit (latest pip-installable version)
- tqdm
- networkx  
- python-tsp
- pandas