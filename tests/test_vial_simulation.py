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


# ── Variability is nucleation-driven: more vials affected over time ────────────
#
# In the corrected model ripening is suppressed, so the per-vial deficit does
# NOT grow with storage. What grows is the FRACTION of vials that have nucleated
# precipitation — this is the "in some samples" / vial-to-vial observation.

def test_mean_deficit_increases_with_storage(stats):
    """Population mean rises because MORE vials nucleate over time."""
    means = [r["mean_deficit_pct"] for r in stats]
    assert all(means[i] <= means[i+1] for i in range(len(means)-1))

def test_fraction_with_deficit_increases_with_storage(stats):
    fracs = [r["frac_with_deficit"] for r in stats]
    assert all(fracs[i] <= fracs[i+1] for i in range(len(fracs)-1))


# ── "In some samples": affected fraction is intermediate, grows with time ──────

def test_fraction_at_1month_is_minority(stats):
    """At 1 month only a minority of vials have nucleated precipitation."""
    row = next(r for r in stats if r["storage_months"] == 1)
    assert row["frac_with_deficit"] < 0.40, (
        f"1-month: {row['frac_with_deficit']:.1%} affected; expected a minority"
    )

def test_fraction_at_6months_is_some_not_all(stats):
    """
    At 6 months a substantial but not universal fraction is affected — this
    brackets the Seeker's 'in some samples' / vial-to-vial observation.
    """
    row = next(r for r in stats if r["storage_months"] == 6)
    f   = row["frac_with_deficit"]
    assert 0.30 <= f <= 0.95, (
        f"6-month: {f:.1%} affected; expected 'some but not all' (0.30–0.95)"
    )

def test_fraction_at_24months_is_large_majority(stats):
    """By 24 months almost all vials have nucleated."""
    row = next(r for r in stats if r["storage_months"] == 24)
    assert row["frac_with_deficit"] > 0.80, (
        f"24-month: {row['frac_with_deficit']:.1%} affected; expected > 80%"
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

def test_deficit_per_nucleated_vial_independent_of_storage():
    """
    With ripening suppressed (corrected pool viscosity), a vial that has
    nucleated shows the SAME per-vial deficit regardless of storage duration —
    the precipitate does not transform/coarsen. (The population effect grows only
    because more vials nucleate; see the fraction tests above.)
    """
    kw = dict(
        local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, storage_T_C=-20.0,
        albumin_g_dL=5.5, cryo_purity=1.0,
        glass_density=0.5, fill_volume_mL=7.5,
    )
    d1 = compute_vial_deficit(storage_months=3,  **kw)
    d2 = compute_vial_deficit(storage_months=12, **kw)
    assert abs(d2 - d1) < 0.005, "per-vial deficit should be ~flat vs storage"

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

def test_nucleated_vial_deficit_matches_module5_band():
    """
    A nucleated vial's per-vial deficit should match Module 5's deterministic
    value (modest, ~1-3% for the representative aggregate size; the magnitude is
    uncertain and experimentally determined — see ACP_AGGREGATE_NM).
    """
    d = compute_vial_deficit(
        storage_months=6, local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, storage_T_C=-20.0,
        albumin_g_dL=5.5, cryo_purity=1.0,
        glass_density=0.5, fill_volume_mL=7.5,
    )
    assert 0.005 <= d <= 0.05, (
        f"nucleated-vial deficit = {d*100:.1f}%; expected modest (0.5–5%)"
    )
