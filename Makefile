PYTHON  := python
PYTEST  := python -m pytest

.PHONY: all setup test figures clean

all: setup test figures
	@echo ""
	@echo "=== make all complete: 139 tests, 11 figures, 6 data tables ==="

# ── Environment ───────────────────────────────────────────────────────────────

setup:
	pip install -e ".[dev]" --quiet
	@echo "✓ Dependencies installed"

# ── Tests ─────────────────────────────────────────────────────────────────────

test:
	$(PYTEST) tests/ -v

# ── Figures ───────────────────────────────────────────────────────────────────

figures:
	@mkdir -p figures data
	$(PYTHON) -m src.freezing_trajectory
	$(PYTHON) -m src.saturation_indices
	$(PYTHON) -m src.supersaturation_map
	$(PYTHON) -m src.ripening_kinetics
	$(PYTHON) -m src.plot_fig04
	$(PYTHON) -m src.plot_fig05
	$(PYTHON) -m src.vial_simulation
	$(PYTHON) -m src.plot_sobol
	$(PYTHON) -m src.plot_fig06
	$(PYTHON) -m src.plot_fig07
	$(PYTHON) -m src.interventions
	$(PYTHON) -m src.plot_fig08_09_10
	$(PYTHON) -m src.phreeqc_runner
	@echo "✓ All figures + data tables regenerated (figures/, data/)"

# ── Cleanup ───────────────────────────────────────────────────────────────────

clean:
	rm -rf figures/*.png figures/*.pdf data/*.csv
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -name "*.pyc" -delete
	@echo "✓ Clean"
