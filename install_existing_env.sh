#!/bin/bash

# Install Molecular Complexity package into existing environment
# This script detects the current environment and installs appropriately

echo "Installing Molecular Complexity package into existing environment..."

# Check if we're in a virtual environment
if [[ -n "$VIRTUAL_ENV" ]]; then
    echo "✓ Detected active virtual environment: $VIRTUAL_ENV"
elif [[ -n "$CONDA_DEFAULT_ENV" ]]; then
    echo "✓ Detected active conda environment: $CONDA_DEFAULT_ENV"
else
    echo "⚠️  No virtual environment detected. Consider activating one first."
    echo "   Examples:"
    echo "   - conda activate your_env"
    echo "   - source your_venv/bin/activate"
    echo ""
    read -p "Continue with system Python? [y/N]: " confirm
    if [[ $confirm != [yY] ]]; then
        echo "Installation cancelled."
        exit 1
    fi
fi

# Check if uv is available
if command -v uv &> /dev/null; then
    echo "Using uv for installation..."
    uv pip install -e .
elif command -v conda &> /dev/null && [[ -n "$CONDA_DEFAULT_ENV" ]]; then
    echo "Using conda for dependencies..."
    conda install -y numpy rdkit pandas
    pip install -e .
elif command -v mamba &> /dev/null && [[ -n "$CONDA_DEFAULT_ENV" ]]; then
    echo "Using mamba for dependencies..."
    mamba install -y numpy rdkit pandas
    pip install -e .
else
    echo "Using pip for installation..."
    pip install -e .
fi

# Test the installation
echo ""
echo "Testing installation..."
if python -c "from molecular_complexity import molecular_complexity; print('✓ Import successful')" 2>/dev/null; then
    echo "✓ Installation successful!"
    echo ""
    echo "Try it out:"
    echo "  python -c \"from molecular_complexity import molecular_complexity; print(molecular_complexity('CCO'))\""
    echo "  molecular-complexity 'CCO' 'CCC'"
    echo "  python examples.py"
else
    echo "✗ Installation failed. Please check error messages above."
    exit 1
fi