"""
Tests for Module 4: supersaturation heatmap grid.
All checks are physical/biological sanity constraints.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.supersaturation_map import (
    compute_si_grid,
    K_GRID,
    PH_GRID,
    T_MAP,
    PHASES_PANEL,
    _comp_at_k,
    _trajectory_k_pH,
)


@pytest.fixture(scope="module")
def grids():
    """Compute the full SI grid once for all tests in this module."""
    return compute_si_grid()


# ── Grid shape and finiteness ─────────────────────────────────────────────────

def test_grid_shape(grids):
    for phase in PHASES_PANEL:
        assert grids[phase].shape == (len(K_GRID), len(PH_GRID))

def test_grid_all_finite(grids):
    for phase in PHASES_PANEL:
        assert np.all(np.isfinite(grids[phase])), f"Non-finite values in {phase} grid"


# ── Monotonicity: higher k → higher SI at fixed pH ───────────────────────────

def test_HAp_SI_increases_with_k(grids):
    """More concentrated → higher SI(HAp) at any fixed pH."""
    ph_idx = int(np.argmin(np.abs(PH_GRID - 7.4)))
    SI_col = grids["HAp"][:, ph_idx]
    assert np.all(np.diff(SI_col) >= -0.01), "SI(HAp) should increase monotonically with k"

def test_ACP_SI_increases_with_k(grids):
    ph_idx = int(np.argmin(np.abs(PH_GRID - 7.4)))
    SI_col = grids["ACP_loose"][:, ph_idx]
    assert np.all(np.diff(SI_col) >= -0.01)


# ── Monotonicity: higher pH → higher SI for CaP phases ───────────────────────

def test_HAp_SI_increases_with_pH(grids):
    """More alkaline → more PO₄³⁻ and OH⁻ → higher SI(HAp)."""
    k_idx = int(np.argmin(np.abs(K_GRID - 5.0)))
    SI_row = grids["HAp"][k_idx, :]
    assert np.all(np.diff(SI_row) >= -0.01)

def test_calcite_SI_increases_with_pH(grids):
    k_idx = int(np.argmin(np.abs(K_GRID - 5.0)))
    SI_row = grids["calcite"][k_idx, :]
    assert np.all(np.diff(SI_row) >= -0.01)


# ── Key reference points ──────────────────────────────────────────────────────

def test_HAp_SI_at_cryo_state(grids):
    """
    At k≈5.58, pH≈7.81, T=−20°C: SI(HAp) should be > physiological (~+5).
    Cryo supersaturation drives the hypothesised precipitation.
    """
    k_idx  = int(np.argmin(np.abs(K_GRID - 5.58)))
    ph_idx = int(np.argmin(np.abs(PH_GRID - 7.81)))
    si = grids["HAp"][k_idx, ph_idx]
    assert si >= 5.0, f"SI(HAp) at cryo state = {si:.2f}; expected ≥ 5"

def test_HAp_supersaturated_over_wide_range(grids):
    """SI(HAp) > 0 for all k > 2 at pH 7.4 — serum is always supersaturated wrt HAp."""
    ph_idx  = int(np.argmin(np.abs(PH_GRID - 7.4)))
    k_mask  = K_GRID > 2.0
    SI_vals = grids["HAp"][k_mask, ph_idx]
    assert np.all(SI_vals > 0), "SI(HAp) should be > 0 for k > 2 at pH 7.4"

def test_brushite_less_supersaturated_than_HAp(grids):
    """Brushite has higher Ksp than HAp → lower SI everywhere."""
    # Check at a representative cryo point
    k_idx  = int(np.argmin(np.abs(K_GRID - 5.58)))
    ph_idx = int(np.argmin(np.abs(PH_GRID - 7.81)))
    assert grids["brushite"][k_idx, ph_idx] < grids["HAp"][k_idx, ph_idx]

def test_ACP_loose_less_supersaturated_than_ACP_tight(grids):
    """ACP_loose has higher Ksp → lower SI than ACP_tight."""
    k_idx  = int(np.argmin(np.abs(K_GRID - 5.58)))
    ph_idx = int(np.argmin(np.abs(PH_GRID - 7.81)))
    assert grids["ACP_loose"][k_idx, ph_idx] < grids["ACP_tight"][k_idx, ph_idx]


# ── Composition scaling ───────────────────────────────────────────────────────

def test_comp_at_k_scales_all_species():
    k = 3.7
    comp = _comp_at_k(k)
    from src.supersaturation_map import BASE_COMP
    for sp, v in BASE_COMP.items():
        assert abs(comp[sp] - v * k) < 1e-10, f"Species {sp} not scaled correctly"


# ── Trajectory sanity ─────────────────────────────────────────────────────────

def test_trajectory_pH_increases_with_co2_loss():
    """At same k, more CO₂ loss → higher pH."""
    k0, pH_0pct  = _trajectory_k_pH(0.0)
    k5, pH_50pct = _trajectory_k_pH(0.5)
    assert np.all(pH_50pct >= pH_0pct - 0.01)

def test_trajectory_k_array_matches_grid():
    k_arr, _ = _trajectory_k_pH(0.0)
    assert len(k_arr) == len(K_GRID)
    assert np.allclose(k_arr, K_GRID)
