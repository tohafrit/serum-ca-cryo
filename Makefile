PYTHON  := .venv/bin/python
PIP     := .venv/bin/pip
PYTEST  := .venv/bin/pytest
JUPYTER := .venv/bin/jupyter

.PHONY: all setup test figures notebooks clean

all: setup test figures

# ── Environment ───────────────────────────────────────────────────────────────

setup: .venv/bin/activate

.venv/bin/activate: pyproject.toml
	python3 -m venv .venv
	$(PIP) install --upgrade pip
	$(PIP) install -e ".[dev]"
	@echo "✓ Environment ready. Activate with: source .venv/bin/activate"

# ── Tests ─────────────────────────────────────────────────────────────────────

test: setup
	$(PYTEST) tests/ -v

# ── Figures ───────────────────────────────────────────────────────────────────

figures: setup
	@mkdir -p figures
	$(PYTHON) src/freezing_trajectory.py
	$(PYTHON) src/saturation_indices.py
	$(PYTHON) src/supersaturation_map.py
	$(PYTHON) src/ripening_kinetics.py
	$(PYTHON) src/vial_simulation.py
	$(PYTHON) src/interventions.py
	@echo "✓ All figures generated in figures/"

# ── Notebooks ────────────────────────────────────────────────────────────────

notebooks: setup
	$(JUPYTER) notebook notebooks/

# ── Cleanup ───────────────────────────────────────────────────────────────────

clean:
	rm -rf figures/*.png figures/*.pdf data/*.csv
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -name "*.pyc" -delete
