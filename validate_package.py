#!/usr/bin/env python3
"""
SARdine Package Validation Script
=================================

This script validates the SARdine package installation and basic functionality.
Run this after cloning the repository to ensure everything is set up correctly.
"""

import os
import sys
import subprocess
import importlib.util
from pathlib import Path

def check_rust_installation():
    """Check if Rust is installed and available."""
    try:
        result = subprocess.run(['cargo', '--version'], 
                              capture_output=True, text=True, check=True)
        print(f"✅ Rust/Cargo: {result.stdout.strip()}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ Rust/Cargo not found. Please install Rust: https://rustup.rs/")
        return False

def check_python_installation():
    """Check Python version and basic requirements."""
    version = sys.version_info
    if version.major != 3 or version.minor < 8:
        print(f"❌ Python 3.8+ required, found {version.major}.{version.minor}")
        return False
    print(f"✅ Python: {version.major}.{version.minor}.{version.micro}")
    return True

def check_project_structure():
    """Validate the project structure."""
    required_dirs = ['src', 'python', 'examples', 'docs', 'tests']
    required_files = ['Cargo.toml', 'pyproject.toml', 'README.md', 'LICENSE']
    
    missing_dirs = [d for d in required_dirs if not os.path.isdir(d)]
    missing_files = [f for f in required_files if not os.path.isfile(f)]
    
    if missing_dirs:
        print(f"❌ Missing directories: {missing_dirs}")
        return False
    if missing_files:
        print(f"❌ Missing files: {missing_files}")
        return False
    
    print("✅ Project structure is valid")
    return True

def check_rust_build():
    """Check if Rust code compiles."""
    try:
        result = subprocess.run(['cargo', 'check'], 
                              capture_output=True, text=True, check=True)
        print("✅ Rust code compiles successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Rust compilation failed: {e.stderr}")
        return False

def check_python_sardine():
    """Check if sardine Python package can be imported."""
    try:
        # Try to import sardine from python/sardine
        spec = importlib.util.spec_from_file_location("sardine", "python/sardine/__init__.py")
        sardine = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sardine)
        print("✅ SARdine Python package can be imported")
        return True
    except Exception as e:
        print(f"❌ SARdine Python import failed: {e}")
        return False

def check_example_files():
    """Check if key example files exist."""
    key_examples = [
        'examples/complete_backscatter_pipeline.py',
        'examples/backscatter_processor.py',
        'examples/complete_calibration_workflow.py'
    ]
    
    missing = [f for f in key_examples if not os.path.isfile(f)]
    if missing:
        print(f"❌ Missing example files: {missing}")
        return False
    
    print("✅ Key example files are present")
    return True

def check_data_directory():
    """Check if data directory exists and has expected structure."""
    data_dir = Path('data')
    if not data_dir.exists():
        print("❌ Data directory not found")
        return False
    
    # Check for at least one .zip file (SLC data)
    zip_files = list(data_dir.glob('*.zip'))
    if not zip_files:
        print("⚠️  No SLC data files found in data/ directory")
        print("   Note: Add Sentinel-1 SLC data to run complete examples")
        return True
    
    print(f"✅ Data directory contains {len(zip_files)} SLC file(s)")
    return True

def main():
    """Run all validation checks."""
    print("🔍 SARdine Package Validation")
    print("=" * 40)
    
    checks = [
        ("Python Installation", check_python_installation),
        ("Rust Installation", check_rust_installation),
        ("Project Structure", check_project_structure),
        ("Rust Compilation", check_rust_build),
        ("Python Package", check_python_sardine),
        ("Example Files", check_example_files),
        ("Data Directory", check_data_directory),
    ]
    
    results = []
    for name, check_func in checks:
        print(f"\n📋 {name}:")
        try:
            success = check_func()
            results.append(success)
        except Exception as e:
            print(f"❌ {name} failed with error: {e}")
            results.append(False)
    
    print("\n" + "=" * 40)
    print("📊 Validation Summary:")
    passed = sum(results)
    total = len(results)
    print(f"✅ {passed}/{total} checks passed")
    
    if passed == total:
        print("\n🎉 All checks passed! SARdine is ready to use.")
        print("\nNext steps:")
        print("1. Add Sentinel-1 SLC data to the data/ directory")
        print("2. Run: python examples/complete_backscatter_pipeline.py")
        print("3. Explore other examples in the examples/ directory")
    else:
        print(f"\n⚠️  {total - passed} checks failed. Please address the issues above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
