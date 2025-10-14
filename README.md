# HamDiv

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

HamDiv provides an implementation of **Hamiltonian diversity**, together with all existing molecular diversity metrics for chemical datasets.

## Installation

### Quick Setup (New Environment)

Run the setup script to create a new virtual environment and install:

```bash
./setup.sh
```

### Install into Existing Environment

#### Using uv (Recommended)

If you have an existing uv environment:
```bash
# Activate your existing environment
source /path/to/your/.venv/bin/activate
# Install the package
uv pip install -e .
```

#### Using conda/mamba

If you have an existing conda or mamba environment:
```bash
# Activate your conda environment
conda activate your_env_name
# or: mamba activate your_env_name

# Install dependencies
conda install numpy rdkit pandas tqdm networkx
# or: mamba install numpy rdkit pandas tqdm networkx

# Install python-tsp via pip (not available in conda)
pip install python-tsp

# Install this package
pip install -e .
```

#### Using regular pip

For any Python environment:
```bash
# Activate your environment
source your_env/bin/activate

# Install this package (dependencies will be installed automatically)
pip install -e .
```

### Manual Installation (New Environment)

1. Create a virtual environment:
```bash
uv venv
source .venv/bin/activate
```

2. Install the package:
```bash
uv pip install -e .
```

### Alternative: Install from requirements

```bash
uv pip install -r requirements.txt
uv pip install -e .
```

## Usage

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
print(f"Richness: {richness}")

int_div = diversity_all(smiles=smiles, mode="IntDiv")
print(f"IntDiv: {int_div}")

circles = diversity_all(smiles=smiles, mode="NCircles-0.7")  # 0.7 is adjustable
print(f"Circles: {circles}")

# Hamiltonian diversity (main contribution)
ham_div = diversity_all(smiles=smiles, mode="HamDiv")
# or directly:
ham_div_direct = HamDiv(smiles=smiles)
print(f"HamDiv: {ham_div}")

# Get all diversity metrics as a vector
diversity_dict = diversity_vector(smiles=smiles)
print("All metrics:", diversity_dict)

# Get all diversity metrics as a numpy array
metric_names, values_array = diversity_vector_array(smiles=smiles)
print("Metric names:", metric_names)
print("Values:", values_array)
```

## Available Diversity Metrics

### Count-based Metrics
- **Richness**: Number of unique molecules
- **FG**: Functional groups diversity  
- **RS**: Ring systems diversity
- **BM**: Bemis-Murcko scaffolds diversity

### Distance-based Metrics
- **IntDiv**: Internal diversity (average pairwise distance)
- **SumDiv**: Sum diversity
- **Diam**: Diameter (maximum pairwise distance)
- **SumDiam**: Sum diameter
- **Bot**: Bottleneck (minimum pairwise distance)
- **SumBot**: Sum bottleneck
- **DPP**: Determinantal point process
- **NCircles-X**: N-circles with threshold X
- **HamDiv**: Hamiltonian diversity (solving TSP on molecular distance graph)

## Key Features

- **Hamiltonian Diversity**: Novel diversity metric based on traveling salesman problem solution
- **Comprehensive Metrics**: Includes all major molecular diversity measures
- **Vector Output**: Get all metrics as a dictionary or numpy array with `diversity_vector()` and `diversity_vector_array()`
- **RDKit Integration**: Uses ECFP fingerprints and Tanimoto similarity
- **TSP Algorithms**: Multiple TSP solving methods (greedy, exact, heuristics)
- **Easy Installation**: Simple setup with `uv` package manager
- **Latest Dependencies**: Uses modern `rdkit` from pip (not rdkit-pypi)

## Dependencies

- Python 3.8+
- numpy
- rdkit (modern pip-installable version)
- tqdm
- networkx
- python-tsp
- pandas

## Additional Features

The `additional_term.py` module provides core functionality for incorporating Hamiltonian diversity into molecular generation scoring functions (Section 4.2 in the paper).

## Example Output

```
=== Individual Diversity Metrics ===
Richness: 5
IntDiv: 0.6396810755759057
Circles: 2.0
HamDiv: 3.0460708271069046

=== All Diversity Metrics as Vector ===
All metrics: {
    'richness': 5,
    'functional_groups': 8,
    'ring_systems': 2,
    'bemis_murcko': 4,
    'internal_diversity': 0.6396810755759057,
    'sum_diversity': 12.793621511518114,
    'diameter': 0.9848484848484849,
    'sum_diameter': 4.090909090909091,
    'bottleneck': 0.3181818181818182,
    'sum_bottleneck': 2.590909090909091,
    'determinantal_point_process': 0.030047067904866817,
    'n_circles': 2.0,
    'hamiltonian_diversity': 3.0460708271069046
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use HamDiv in your research, please cite the original paper describing Hamiltonian diversity.
