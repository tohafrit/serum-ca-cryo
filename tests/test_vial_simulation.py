"""
Tests for Module 6: Monte Carlo vial-to-vial simulation.
All checks are physical/biological sanity constraints.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.vial_simulation import (
    run_monte_carlo,
    vial_statistics,
    compute_vial_deficit,
    STORAGE_MONTHS,
    THRESHOLD,
    BASE_K,
    NUCLEATION_MEDIAN_DAYS,
    N_VIALS,
)


@pytest.fixture(scope="module")
def mc_result():
    return run_monte_carlo(seed=42)


@pytest.fixture(scope="module")
def stats(mc_result):
    return vial_statistics(mc_result["deficits"])


# ── Output shape and finiteness ────────────────────────────────────────────────

def test_deficits_shape(mc_result):
    assert mc_result["deficits"].shape == (N_VIALS, len(STORAGE_MONTHS))

def test_deficits_all_finite(mc_result):
    assert np.all(np.isfinite(mc_result["deficits"]))

def test_deficits_non_negative(mc_result):
    assert np.all(mc_result["deficits"] >= 0)

def test_deficits_bounded(mc_result):
    assert np.all(mc_result["deficits"] <= 1.0)


# ── Monotonicity: deficit grows with storage time ─────────────────────────────

def test_mean_deficit_increases_with_storage(stats):
    means = [r["mean_deficit_pct"] for r in stats]
    assert all(means[i] <= means[i+1] for i in range(len(means)-1))

def test_fraction_above_threshold_increases_with_storage(stats):
    fracs = [r["fraction_above_4pct"] for r in stats]
    assert all(fracs[i] <= fracs[i+1] for i in range(len(fracs)-1))


# ── Critical threshold tests (Seeker's observations) ─────────────────────────

def test_fraction_above_4pct_at_1month_below_5pct(stats):
    """At 1 month, almost no vials should exceed 4% deficit."""
    row = next(r for r in stats if r["storage_months"] == 1)
    assert row["fraction_above_4pct"] < 0.05, (
        f"1-month: {row['fraction_above_4pct']:.1%} above 4%; expected < 5%"
    )

def test_fraction_above_4pct_at_6months_in_range(stats):
    """
    CRITICAL: At 6 months, 25–75% of vials should exceed 4%.
    Seeker says 'in some samples' → must not be 0% or 100%.
    This brackets the Seeker's vial-to-vial observation.
    """
    row = next(r for r in stats if r["storage_months"] == 6)
    f   = row["fraction_above_4pct"]
    assert 0.25 <= f <= 0.75, (
        f"6-month: {f:.1%} vials above 4%; expected 25–75% "
        "(Seeker: 'vial-to-vial dependent', 'in some samples')"
    )

def test_fraction_above_4pct_at_24months_above_85pct(stats):
    """At 24 months, the large majority of vials should be affected."""
    row = next(r for r in stats if r["storage_months"] == 24)
    assert row["fraction_above_4pct"] > 0.85, (
        f"24-month: {row['fraction_above_4pct']:.1%} above 4%; expected > 85%"
    )


# ── Sobol sanity (run separately; test CSV existence + S1 sum) ────────────────

def test_sobol_csv_exists():
    csv_path = Path(__file__).parent.parent / "data" / "module6_sobol_indices.csv"
    assert csv_path.exists(), "Run plot_sobol.py to generate module6_sobol_indices.csv"

def test_sobol_nucleation_delay_dominates():
    """Nucleation delay must be the top-ranked parameter by total-order ST."""
    import csv
    csv_path = Path(__file__).parent.parent / "data" / "module6_sobol_indices.csv"
    if not csv_path.exists():
        pytest.skip("Sobol CSV not yet generated")
    with open(csv_path) as f:
        rows = list(csv.DictReader(f))
    # Rows are sorted by ST descending
    top_param = rows[0]["parameter"]
    assert top_param == "log10_nuc_delay_days", (
        f"Top Sobol parameter is {top_param}; expected nucleation delay"
    )


# ── compute_vial_deficit sanity ───────────────────────────────────────────────

def test_deficit_zero_when_nucleation_not_yet_occurred():
    """If nucleation delay > storage time, deficit must be 0."""
    d = compute_vial_deficit(
        storage_months        = 1.0,
        local_k_factor        = 1.0,
        nucleation_delay_days = 999.0,   # much longer than 30 days
        particle_size_nm      = 50.0,
        storage_T_C           = -20.0,
        albumin_g_dL          = 5.5,
        cryo_purity           = 1.0,
        glass_density         = 0.5,
        fill_volume_mL        = 7.5,
    )
    assert d == 0.0

def test_deficit_increases_with_effective_storage_time():
    """Longer storage (same delay) → more HAp → higher deficit."""
    kw = dict(
        local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, storage_T_C=-20.0,
        albumin_g_dL=5.5, cryo_purity=1.0,
        glass_density=0.5, fill_volume_mL=7.5,
    )
    d1 = compute_vial_deficit(storage_months=3,  **kw)
    d2 = compute_vial_deficit(storage_months=12, **kw)
    assert d2 > d1

def test_deficit_increases_with_warmer_storage_T():
    """
    Warmer storage (at fixed k) → faster Arrhenius kinetics → more ACP→HAp
    transformation → higher 60-min deficit.  The kinetic effect dominates
    when cryoconcentration factor k is held constant.
    """
    kw = dict(
        storage_months=12, local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, albumin_g_dL=5.5, cryo_purity=1.0,
        glass_density=0.5, fill_volume_mL=7.5,
    )
    d_cold = compute_vial_deficit(storage_T_C=-22.0, **kw)
    d_warm = compute_vial_deficit(storage_T_C=-18.0, **kw)
    assert d_warm >= d_cold   # warmer = faster Arrhenius = more HAp at fixed k

def test_median_nucleated_vials_close_to_module5():
    """
    Vials with zero nucleation delay should give deficit close to Module 5's
    deterministic prediction (~5.4% at 6 months).  Checks model consistency.
    """
    d = compute_vial_deficit(
        storage_months=6, local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, storage_T_C=-20.0,
        albumin_g_dL=5.5, cryo_purity=1.0,
        glass_density=0.5, fill_volume_mL=7.5,
    )
    assert 0.03 <= d <= 0.10, (
        f"Module 5-equivalent deficit = {d*100:.1f}%; expected 3–10%"
    )
