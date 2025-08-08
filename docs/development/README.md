# SARdine - SAR Processing Pipeline

⚠️ **IMPORTANT DISCLAIMER: THIS PACKAGE IS NOT PRODUCTION READY** ⚠️

## Development Status

🚧 **EXPERIMENTAL / DEVELOPMENT PHASE** 🚧

This package is currently in **early development** and should be considered a **playground/experimental tool only**. 

**DO NOT USE IN PRODUCTION ENVIRONMENTS**

### Current Status:
- ✅ Core SAR processing pipeline implemented
- ✅ Optimized terrain correction integrated  
- ✅ 14-step backscatter processing workflow
- ✅ CLI interface functional
- ❌ **NOT** production ready
- ❌ **NOT** fully tested
- ❌ **NOT** validated for scientific accuracy
- ❌ **NOT** suitable for operational use

### Known Limitations:
- Limited error handling in edge cases
- Incomplete test coverage
- Performance not optimized for large datasets
- Scientific validation incomplete
- Documentation incomplete
- No stability guarantees

## What is SARdine?

SARdine is an experimental SAR (Synthetic Aperture Radar) processing pipeline built with Rust and Python, designed to process Sentinel-1 data. It provides a complete 14-step workflow for SAR backscatter processing.

### Key Features (Experimental):
- **14-step SAR processing pipeline**
- **Optimized terrain correction** (10-20x performance improvement)
- **Real data processing** (no synthetic fallbacks)
- **Scientific validation enforcement**
- **CLI interface**
- **Rust/Python integration**

## Installation (Development Only)

```bash
# Clone the repository
git clone https://github.com/SteveMHill/SARdine.git
cd SARdine

# This package is in development - installation process may change
# Currently requires manual setup
```

## Usage (Experimental)

⚠️ **WARNING: Use only for testing and development purposes** ⚠️

```bash
# Example command (subject to change)
python3 -m sardine.cli backscatter input.zip output_dir --polarization VV --resolution 10
```

## Project Structure

```
SARdine/
├── SARdine/           # Main package
│   ├── python/        # Python modules
│   └── rust/          # Rust implementations
├── docs/              # Documentation (incomplete)
├── scripts/           # Development scripts
├── tests/             # Test files (incomplete coverage)
└── README.md          # This file
```

## Development

This is a **research/development project**. Contributions are welcome but understand that:

- APIs may change without notice
- No backward compatibility guarantees
- Code quality varies across modules
- Performance characteristics not finalized

## Disclaimer

**THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND**

- Use at your own risk
- No guarantee of scientific accuracy
- May produce incorrect results
- Not suitable for any operational use
- Development quality code only

## Future Plans

- Complete scientific validation
- Comprehensive testing
- Performance optimization
- Production-ready releases
- Proper documentation
- Stability guarantees

## Contact

This is an experimental project. For questions about the development:
- Check the issues on GitHub
- Understand this is not production software

---

**Remember: This is experimental software - NOT FOR PRODUCTION USE**
