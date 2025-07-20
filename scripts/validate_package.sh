#!/bin/bash

# =============================================================================
# SARdine Package Validation Script
# =============================================================================
# Validates that SARdine is ready for GitHub upload

set -e

echo "🔍 Validating SARdine Package for GitHub Upload..."
echo "=================================================="

# Check essential files exist
echo "📋 Checking essential files..."
essential_files=(
    "README.md"
    "LICENSE" 
    "Cargo.toml"
    "pyproject.toml"
    "build.sh"
    "logo.png"
    "CONTRIBUTING.md"
    ".gitignore"
)

for file in "${essential_files[@]}"; do
    if [[ -f "$file" ]]; then
        echo "  ✅ $file"
    else
        echo "  ❌ $file (missing)"
        exit 1
    fi
done

# Check directory structure
echo ""
echo "📁 Checking directory structure..."
essential_dirs=(
    "src"
    "python"
    "docs"
    "examples"
    "tests"
    "development"
)

for dir in "${essential_dirs[@]}"; do
    if [[ -d "$dir" ]]; then
        echo "  ✅ $dir/"
    else
        echo "  ❌ $dir/ (missing)"
        exit 1
    fi
done

# Check that development files are moved
echo ""
echo "🧹 Checking development files are organized..."
if [[ -f "test_*.py" ]] || [[ -f "validate_*.py" ]] || [[ -f "*_analysis.py" ]]; then
    echo "  ⚠️  Some development files may still be in root directory"
    echo "  Consider moving them to development/ directory"
else
    echo "  ✅ Development files properly organized"
fi

# Check build capability
echo ""
echo "🔨 Testing build capability..."
if command -v cargo &> /dev/null; then
    echo "  ✅ Rust/Cargo available"
    
    if command -v python3 &> /dev/null || command -v python &> /dev/null; then
        echo "  ✅ Python available"
        
        if command -v maturin &> /dev/null; then
            echo "  ✅ Maturin available"
        else
            echo "  ⚠️  Maturin not installed (will be installed during build)"
        fi
    else
        echo "  ❌ Python not available"
        exit 1
    fi
else
    echo "  ❌ Rust/Cargo not available"
    exit 1
fi

# Check git status
echo ""
echo "📦 Checking git status..."
if [[ -d ".git" ]]; then
    echo "  ✅ Git repository initialized"
    
    # Check for large files that shouldn't be committed
    echo "  🔍 Checking for large files..."
    large_files=$(find . -type f -size +10M -not -path "./.git/*" -not -path "./target/*" -not -path "./.venv/*" 2>/dev/null || true)
    if [[ -n "$large_files" ]]; then
        echo "  ⚠️  Large files found (consider adding to .gitignore):"
        echo "$large_files" | sed 's/^/    /'
    else
        echo "  ✅ No large files found"
    fi
else
    echo "  ⚠️  Not a git repository"
fi

echo ""
echo "✅ Package validation completed!"
echo "=================================================="
echo ""
echo "🚀 Ready for GitHub upload! Next steps:"
echo "  1. git add ."
echo "  2. git commit -m 'chore: prepare package for public release'"
echo "  3. git push origin main"
echo "  4. Create release on GitHub"
echo ""
echo "📚 Don't forget to:"
echo "  • Add GitHub repository description"
echo "  • Add topics/tags (rust, python, sar, remote-sensing, sentinel-1)"
echo "  • Enable GitHub Pages for documentation (if desired)"
echo "  • Add GitHub Actions for CI/CD (optional)"
