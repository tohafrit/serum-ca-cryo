"""
Trivial smoke tests that verify the project scaffold and Python environment
are correctly set up before any scientific modules are loaded.
"""

import sys
import importlib


def test_python_version():
    """Require Python 3.11+ for match-statement and tomllib support."""
    assert sys.version_info >= (3, 11), (
        f"Python 3.11+ required, got {sys.version_info.major}.{sys.version_info.minor}"
    )


def test_core_imports():
    """All mandatory scientific libraries must be importable."""
    for lib in ("numpy", "scipy", "matplotlib", "pandas"):
        mod = importlib.import_module(lib)
        assert mod is not None, f"Failed to import {lib}"


def test_salib_import():
    """SALib is needed for Sobol sensitivity analysis in Module 6."""
    salib = importlib.import_module("SALib")
    assert salib is not None


def test_figures_dir_exists():
    """figures/ directory must exist so matplotlib can save output there."""
    from pathlib import Path
    figures = Path(__file__).parent.parent / "figures"
    assert figures.is_dir(), "figures/ directory missing — run 'make setup' or create it manually"


def test_data_dir_exists():
    """data/ directory must exist for CSV exports."""
    from pathlib import Path
    data = Path(__file__).parent.parent / "data"
    assert data.is_dir(), "data/ directory missing"


def test_basic_arithmetic():
    """Sanity check that numpy is functional."""
    import numpy as np
    result = np.dot([3.0, 4.0], [3.0, 4.0])
    assert abs(result - 25.0) < 1e-12
