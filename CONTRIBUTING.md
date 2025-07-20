# Contributing to SARdine

We welcome contributions to SARdine! This document provides guidelines for contributing to the project.

## 🚨 Alpha Release Status

SARdine is currently in alpha development. APIs may change, and some features are still being implemented. Please keep this in mind when contributing.

## 🔄 Development Workflow

### 1. Setting Up Development Environment

```bash
# Fork the repository on GitHub
git clone https://github.com/your-username/SARdine.git
cd SARdine

# Set up development environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install maturin pytest
maturin develop

# Verify setup
python -c "import sardine; print('Setup successful!')"
```

### 2. Making Changes

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**:
   - Follow existing code style
   - Add tests for new functionality
   - Update documentation as needed

3. **Test your changes**:
   ```bash
   # Run Rust tests
   cargo test
   
   # Run Python tests
   python -m pytest tests/
   
   # Test the build
   ./build.sh
   ```

### 3. Submitting Changes

1. **Commit your changes**:
   ```bash
   git add .
   git commit -m "feat: describe your changes"
   ```

2. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

3. **Create a Pull Request** on GitHub

## 📝 Code Style Guidelines

### Rust Code
- Follow `rustfmt` formatting (run `cargo fmt`)
- Use `clippy` for linting (run `cargo clippy`)
- Add documentation for public functions
- Include unit tests for new functionality

### Python Code
- Follow PEP 8 style guidelines
- Use type hints where possible
- Add docstrings for public functions
- Include tests in the `tests/` directory

## 🧪 Testing

### Running Tests

```bash
# Rust unit tests
cargo test

# Python tests
python -m pytest tests/ -v

# Integration tests
python examples/basic_workflow.py
```

### Adding Tests

- **Rust tests**: Add `#[cfg(test)]` modules in relevant source files
- **Python tests**: Add test files in `tests/` directory following pytest conventions
- **Integration tests**: Add examples in `examples/` directory

## 📚 Documentation

### Code Documentation
- Rust: Use `///` for public function documentation
- Python: Use docstrings following Google/NumPy style

### User Documentation
- Update README.md for user-facing changes
- Add examples to `examples/` directory
- Update API documentation in `docs/`

## 🐛 Reporting Issues

### Bug Reports
Please include:
- SARdine version
- Operating system and Python version
- Minimal code example to reproduce the issue
- Expected vs actual behavior
- Error messages and stack traces

### Feature Requests
Please include:
- Description of the feature
- Use case and motivation
- Proposed API (if applicable)
- Willingness to implement

## 📋 Development Priorities

### Current Focus Areas
1. **Core Algorithm Stability**: Improving reliability of SAR processing algorithms
2. **API Stabilization**: Finalizing Python and CLI interfaces before beta
3. **Performance Optimization**: Rust backend optimizations
4. **Documentation**: Comprehensive user and developer documentation
5. **Testing**: Expanding test coverage

### Future Roadmap
- Polarimetric processing capabilities
- Time series analysis tools
- Advanced terrain correction algorithms
- Cloud processing integration

## 🤝 Community Guidelines

- Be respectful and inclusive
- Help others learn and contribute
- Share knowledge and best practices
- Report issues constructively
- Follow the [Code of Conduct](CODE_OF_CONDUCT.md)

## 📞 Getting Help

- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: Questions and general discussion
- **Documentation**: Check `docs/` directory and README

## 📄 License

By contributing to SARdine, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to SARdine! 🌟
