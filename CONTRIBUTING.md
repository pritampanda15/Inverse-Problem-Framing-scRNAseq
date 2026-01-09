# Contributing to InverseSC

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/yourusername/inverse-problem-scrna.git
   cd inverse-problem-scrna
   ```
3. Create a development environment:
   ```bash
   conda env create -f environment.yml
   conda activate inverse-sc
   pip install -e ".[dev]"
   ```

## Development Workflow

1. Create a new branch for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes

3. Run tests:
   ```bash
   pytest tests/
   ```

4. Format code:
   ```bash
   black inverse_sc/
   flake8 inverse_sc/
   ```

5. Commit and push:
   ```bash
   git add .
   git commit -m "Description of changes"
   git push origin feature/your-feature-name
   ```

6. Open a Pull Request

## Code Style

- Follow PEP 8 guidelines
- Use type hints where appropriate
- Write docstrings for all public functions (Google style)
- Keep functions focused and modular

## Testing

- Write tests for new features
- Ensure all tests pass before submitting PR
- Aim for >80% code coverage

## Documentation

- Update docstrings for modified functions
- Add examples to notebooks if relevant
- Update API reference if adding new public functions

## Areas for Contribution

### High Priority
- [ ] Batch effect modeling in measurement operator
- [ ] Multi-modal integration (CITE-seq, spatial)
- [ ] Improved scalability for >100k cells
- [ ] Real data benchmarks

### Medium Priority
- [ ] Additional biological priors (GRN, pathways)
- [ ] Improved visualization tools
- [ ] Integration with other frameworks (scVI, Pegasus)

### Low Priority
- [ ] Alternative inference methods (HMC, importance sampling)
- [ ] Web interface for interactive analysis
- [ ] Cloud deployment tools

## Questions?

Open an issue or contact: pritam@stanford.edu

## Code of Conduct

Be respectful and inclusive. We welcome contributors from all backgrounds.
