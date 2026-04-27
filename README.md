# PolyQ

PolyQ is a Python library that provides a novel approach to simulate Quantum Circuits using Boolean Polynomials.

## Features
- Simulate quantum circuits efficiently.
- Built-in support for Boolean polynomial operations.
- Compatible with Python 3.11 and above.

## Installation

To install PolyQ, clone the repository and build the library:

```bash
# Clone the repository
git clone https://github.com/QSDAL-IITR/PolyQ.git
cd PolyQ

pip install -e .


# Install dependencies
pip install -r requirements.txt

# Build the library
python -m build

# Install the library
pip install dist/PolyQ-0.9.0-py3-none-any.whl
```

## Usage

Import the library and use its modules:

```python
from PolyQ import engine, branching

# Example usage
engine.run_simulation()
branching.perform_branching()
```

## Testing

To run the tests, use:

```bash
pytest PolyQ/test_engine
```

## Clearing Cache

If you encounter issues related to Python's cached files, you can clear them by removing all `__pycache__` directories. Run the following command in the root directory of the project:

```bash
find . -name "__pycache__" -exec rm -r {} +
```

## License

This project is licensed under the terms of the license specified in the `LICENSE` file.