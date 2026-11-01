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
```

via make 

```bash
make build
```

or manually

```bash
#create venv
python3 -m venv polyq-dev
source polyq-dev/bin/activate
pip install -e .
# Install dependencies
pip install -r requirements.txt
# Build the library
python -m build
# Install the library
pip install dist/PolyQ-0.2.0-py3-none-any.whl 
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
pytest
```

## Clean up

If you want clean up after the work is done

```bash
deactivate
rm -rf polyq-dev
rm -rf build/ dist/ .eggs/ __pycache__/ .pytest_cache/ __pycache__
find . -type d -name "__pycache__" -exec rm -rf {} +
rm -rf *.egg-info
```

or via make
```bash
make clean
```

## License

This project is licensed under the terms of the license specified in the `LICENSE` file.