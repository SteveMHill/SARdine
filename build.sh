#!/bin/bash

# SARdine Build Script
# Builds the Rust backend and Python package

set -e

echo "🔨 Building SARdine..."

# Check if maturin is installed
if ! command -v maturin &> /dev/null; then
    echo "❌ maturin not found. Installing..."
    pip install maturin
fi

# Clean previous builds
echo "🧹 Cleaning previous builds..."
rm -rf target/
rm -rf dist/
find . -name "*.so" -delete
find . -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true

# Build the Rust extension
echo "🦀 Building Rust extension..."
maturin develop --release

# Run Python tests
echo "🐍 Testing Python interface..."
python -c "
import sardine
print(f'✅ SARdine {sardine.__version__} loaded successfully')
print(f'📦 Available modules: {list(sardine.__dict__.keys())}')
"

echo "✅ Build completed successfully!"
echo ""
echo "Try the CLI:"
echo "  sardine --help"
echo "  sardine info /path/to/sentinel1.zip"
