"""Simple test to demonstrate the diversity vector functionality."""

import warnings
warnings.filterwarnings('ignore')  # Suppress RDKit deprecation warnings for cleaner output

from hamdiv import diversity_vector, diversity_vector_array

# Simple test molecules
smiles = ['CCO', 'CCC', 'CCCC', 'c1ccccc1']

print("=== Testing diversity_vector() ===")
diversity_dict = diversity_vector(smiles=smiles, disable=True)
print("All diversity metrics as dictionary:")
for metric, value in diversity_dict.items():
    print(f"  {metric}: {value:.4f}")

print("\n=== Testing diversity_vector_array() ===")
metric_names, values_array = diversity_vector_array(smiles=smiles, disable=True)
print(f"Number of metrics: {len(metric_names)}")
print(f"Array shape: {values_array.shape}")
print(f"Array dtype: {values_array.dtype}")

print("\nFirst 5 metrics:")
for i in range(5):
    print(f"  {metric_names[i]}: {values_array[i]:.4f}")

print(f"\nHamiltonian diversity: {values_array[-1]:.4f}")