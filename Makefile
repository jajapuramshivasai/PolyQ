# Makefile for PolyQ

.PHONY: build test clean all

# Variables
VENV_DIR = polyq-dev
PYTHON = python3
PIP = pip

# Combined build target: creates venv, installs deps, and builds the library
build:
	@echo "Creating virtual environment..."
	$(PYTHON) -m venv $(VENV_DIR)
	@echo "Installing dependencies..."
	$(VENV_DIR)/bin/$(PIP) install -e .
	$(VENV_DIR)/bin/$(PIP) install -r requirements.txt
	@echo "Building the library..."
	$(VENV_DIR)/bin/$(PYTHON) -m build
	@echo "Build complete. Check the dist/ directory for the .whl file."
	@echo "To use this environment later, run: source $(VENV_DIR)/bin/activate"

# Runs the test suite across all files using the venv's pytest
test:
	$(VENV_DIR)/bin/pytest

# Cleans up all generated files, caches, and the virtual environment
clean:
	rm -rf $(VENV_DIR)
	rm -rf build/ dist/ .eggs/ .pytest_cache/
	find . -type d -name "__pycache__" -exec rm -rf {} +
	rm -rf *.egg-info
	@echo "Cleanup complete."

# Runs the full pipeline: clean, build (which includes venv and install), and test
all: clean build test