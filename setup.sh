#!/bin/bash

# HamDiv Installation Script
# This script sets up the HamDiv package for installation with uv

echo "Setting up HamDiv package..."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "Error: uv is not installed. Please install uv first:"
    echo "  curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# Create virtual environment with latest Python (compatible with new rdkit)
echo "Creating virtual environment..."
uv venv

# Activate virtual environment and install package
echo "Installing HamDiv package..."
source .venv/bin/activate
uv pip install -e .

echo ""
echo "Installation complete! To use HamDiv:"
echo "1. Activate the environment: source .venv/bin/activate"
echo "2. Import in Python: from hamdiv import diversity_all, HamDiv"
echo ""
echo "Test the installation by running: python example.py"