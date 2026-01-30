#!/bin/bash
# Install PyYAML dependency for the workflow

echo "Installing PyYAML for ASH xTB workflow..."

# Try to detect if we're in a conda environment
if [ -n "$CONDA_DEFAULT_ENV" ]; then
    echo "Detected conda environment: $CONDA_DEFAULT_ENV"
    echo "Installing via conda..."
    conda install -y pyyaml
elif command -v conda &> /dev/null; then
    echo "Conda is available. Installing PyYAML in 'ash' environment..."
    conda install -n ash -y pyyaml
else
    echo "Conda not detected. Installing via pip..."
    pip install pyyaml
fi

# Verify installation
python3 -c "import yaml; print('✓ PyYAML installed successfully! Version:', yaml.__version__)" 2>&1

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Installation complete!"
    echo "You can now run: python ash_xtb_explicit_workflow_smoketest.py"
else
    echo ""
    echo "✗ Installation failed. Please install manually:"
    echo "  conda install pyyaml"
    echo "  or"
    echo "  pip install pyyaml"
    exit 1
fi
