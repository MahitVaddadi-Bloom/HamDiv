"""Example usage of HamDiv package."""

from hamdiv import diversity_all, HamDiv, diversity_vector, diversity_vector_array

# Example SMILES from the README
smiles = [
    'Cc1cc(C(O)CNC2(C)CC2S)ccc1O',
    'CCc1cc(C(O)CNC(C)(C)C)ccc1O',
    'CCCC(C)NCC(O)c1ccc(O)c(CO)c1',
    'CNCC(S)c1ccc(O)c(CO)c1',
    'CNCC(C)(C)C(O)c1ccc(S)c(CO)c1'
]

print("=== Individual Diversity Metrics ===")
# Calculate different diversity metrics
richness = diversity_all(smiles=smiles, mode="Richness")
print(f"Richness: {richness}")

int_div = diversity_all(smiles=smiles, mode="IntDiv")
print(f"IntDiv: {int_div}")

circles = diversity_all(smiles=smiles, mode="NCircles-0.7")  # 0.7 is an adjustable hyper-parameter
print(f"Circles: {circles}")

ham_div = diversity_all(smiles=smiles, mode="HamDiv")
ham_div_direct = HamDiv(smiles=smiles)
print(f"HamDiv: {ham_div}, HamDiv (direct): {ham_div_direct}")

print("\n=== All Diversity Metrics as Vector ===")
# Get all metrics as a dictionary
diversity_dict = diversity_vector(smiles=smiles)
print("Diversity metrics dictionary:")
for metric, value in diversity_dict.items():
    print(f"  {metric}: {value}")

print("\n=== All Diversity Metrics as NumPy Array ===")
# Get all metrics as a numpy array
metric_names, values_array = diversity_vector_array(smiles=smiles)
print("Metric names:", metric_names)
print("Values array:", values_array)
print("Array shape:", values_array.shape)