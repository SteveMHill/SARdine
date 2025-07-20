#!/bin/bash

# =============================================================================
# SARdine Build Script
# =============================================================================
# Builds the Rust backend and Python package for SARdine
# Supports both development and production builds

set -e

echo "🔨 Building SARdine SAR Processing Library..."
echo "================================================"

# Check Python availability
if ! command -v python3 &> /dev/null && ! command -v python &> /dev/null; then
    echo "❌ Python not found. Please install Python 3.8+ and try again."
    exit 1
fi

# Use python3 if available, fallback to python
PYTHON_CMD="python3"
if ! command -v python3 &> /dev/null; then
    PYTHON_CMD="python"
fi

echo "🐍 Using Python: $($PYTHON_CMD --version)"

# Check Rust availability
if ! command -v cargo &> /dev/null; then
    echo "❌ Rust/Cargo not found. Please install Rust from https://rustup.rs/"
    exit 1
fi

echo "🦀 Using Rust: $(cargo --version)"

# Check if maturin is installed
if ! command -v maturin &> /dev/null; then
    echo "📦 Installing maturin (Rust-Python bridge)..."
    $PYTHON_CMD -m pip install maturin
fi

echo "🔧 Maturin version: $(maturin --version)"

# Clean previous builds
echo ""
echo "🧹 Cleaning previous builds..."
rm -rf target/wheels/ target/release/ 2>/dev/null || true
rm -rf dist/ build/ *.egg-info/ 2>/dev/null || true
find . -name "*.so" -delete 2>/dev/null || true
find . -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true

# Build the Rust extension
echo ""
echo "🦀 Building Rust extension with optimizations..."
maturin develop --release

# Verify Python interface
echo ""
echo "🐍 Testing Python interface..."
$PYTHON_CMD -c "
import sardine
print(f'✅ SARdine {sardine.__version__} loaded successfully')
print(f'📦 Available modules: {sorted([k for k in sardine.__dict__.keys() if not k.startswith(\"_\")])}')
"

echo ""
echo "✅ Build completed successfully!"
echo "================================================"
echo ""
echo "🚀 Quick Start:"
echo "  # Python API"
echo "  python3 -c \"import sardine; print('SARdine ready!')\""
echo ""
echo "  # CLI Tools"
echo "  sardine --help"
echo "  sardine info /path/to/sentinel1.zip"
echo ""
echo "📚 Documentation: docs/"
echo "🔍 Examples: examples/"
echo "🧪 Tests: python -m pytest tests/"
