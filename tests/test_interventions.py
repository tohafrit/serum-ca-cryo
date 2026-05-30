"""
Tests for Module 7: intervention scenarios.
All tests verify physical sanity and proposal-level claims.

Baseline = sealed vial (pH 7.81). `loose_seal` (pH 8.0) is a RISK state, not an
intervention, so it is excluded from "every intervention improves" checks.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
import numpy as np
from src.interventions import (
    run_all_scenarios, scenario_stats,
    compute_vial_deficit_scenario,
    SCENARIOS, INTERVENTION_NAMES, STORAGE_MONTHS, THRESHOLD,
)


@pytest.fixture(scope="module")
def all_results():
    return run_all_scenarios(seed=42)


@pytest.fixture(scope="module")
def stats(all_results):
    return {r["scenario"]: r for r in scenario_stats(all_results)}


# ── Interventions reduce deficit vs baseline ──────────────────────────────────

def test_each_intervention_reduces_mean_deficit_at_12mo(stats):
    """Every genuine intervention must lower the 12-month mean deficit."""
    baseline = stats["baseline"]["mean_deficit_12mo_pct"]
    for name in INTERVENTION_NAMES:
        val = stats[name]["mean_deficit_12mo_pct"]
        assert val <= baseline, (
            f"{name}: mean deficit {val:.1f}% > baseline {baseline:.1f}%"
        )

def test_each_intervention_reduces_fraction_at_6mo(stats):
    """Every genuine intervention must reduce fraction >4% at 6 months."""
    baseline = stats["baseline"]["frac_with_deficit_6mo"]
    for name in INTERVENTION_NAMES:
        val = stats[name]["frac_with_deficit_6mo"]
        assert val <= baseline, (
            f"{name}: frac>4% at 6mo = {val:.3f} > baseline {baseline:.3f}"
        )


# ── loose-seal is a risk state: it must be WORSE than the sealed baseline ─────

def test_loose_seal_worse_than_baseline(stats):
    """
    CO₂ outgassing through a leaking closure raises pH (7.81→8.0) and must
    make the deficit worse — this is the rationale for tight sealing/degassing.
    """
    f_base  = stats["baseline"]["frac_with_deficit_6mo"]
    f_loose = stats["loose_seal"]["frac_with_deficit_6mo"]
    assert f_loose >= f_base, (
        f"loose_seal {f_loose:.3f} should be ≥ sealed baseline {f_base:.3f}"
    )


# ── Vortex is the most direct single intervention at 6 months ─────────────────

def test_vortex_meaningfully_reduces_6mo_fraction(stats):
    """
    30-s vortex attacks dissolution directly (Noyes-Whitney) and should cut the
    6-month affected fraction by at least 20% relative to baseline.
    """
    f_base   = stats["baseline"]["frac_with_deficit_6mo"]
    f_vortex = stats["+vortex_30s"]["frac_with_deficit_6mo"]
    assert f_vortex <= 0.8 * f_base, (
        f"vortex_30s {f_vortex:.3f} should be ≤ 80% of baseline {f_base:.3f}"
    )


# ── CRF helps but does not eliminate the deficit on its own ───────────────────

def test_crf_reduces_but_insufficient_alone_at_12mo(stats):
    """CRF lowers the 12-month fraction vs baseline but cannot eliminate it."""
    baseline = stats["baseline"]["frac_with_deficit_12mo"]
    crf      = stats["+crf_2C"]["frac_with_deficit_12mo"]
    assert crf < baseline, "CRF should reduce 12-mo fraction vs baseline"
    assert crf > 0.05,     "CRF alone should not fully eliminate the 12-mo deficit"


# ── Combined intervention headline claims ─────────────────────────────────────

def test_combined_6mo_below_5pct(stats):
    """Combined intervention reduces 6-month fraction >4% to <5%."""
    f = stats["+combined"]["frac_with_deficit_6mo"]
    assert f < 0.05, (
        f"Combined: {f:.1%} vials >4% at 6 months; expected <5%"
    )

def test_combined_reduces_12mo_substantially(stats):
    """Combined intervention reduces 12-month fraction >4% by ≥30 pp absolute."""
    baseline = stats["baseline"]["frac_with_deficit_12mo"]
    combined = stats["+combined"]["frac_with_deficit_12mo"]
    reduction_pp = baseline - combined
    assert reduction_pp >= 0.30, (
        f"Combined reduces 12-mo fraction by {reduction_pp:.2f} pp; expected ≥0.30"
    )


# ── Formulation neutrality (sanity) ──────────────────────────────────────────

def test_formulation_chemistry_unchanged_by_vortex():
    """Vortex changes only thaw_h; chemical/process composition is unchanged."""
    baseline_sc = SCENARIOS["baseline"]
    vortex_sc   = SCENARIOS["+vortex_30s"]
    assert vortex_sc["ph"]       == baseline_sc["ph"]
    assert vortex_sc["k_sig"]    == baseline_sc["k_sig"]
    assert vortex_sc["nuc_mult"] == baseline_sc["nuc_mult"]
    assert vortex_sc["thaw_h"] < baseline_sc["thaw_h"]

def test_formulation_chemistry_unchanged_by_crf():
    """CRF changes only k_sig and nuc_mult; pH and thaw are unchanged."""
    baseline_sc = SCENARIOS["baseline"]
    crf_sc      = SCENARIOS["+crf_2C"]
    assert crf_sc["ph"]    == baseline_sc["ph"]
    assert crf_sc["thaw_h"] == baseline_sc["thaw_h"]
    assert crf_sc["k_sig"]  < baseline_sc["k_sig"]


# ── compute_vial_deficit_scenario sanity ──────────────────────────────────────

def test_per_vial_deficit_independent_of_ph():
    """
    In the corrected model pH acts on the AFFECTED FRACTION (via nucleation),
    not on a nucleated vial's per-vial magnitude. A forced-nucleated vial shows
    the same deficit regardless of pH. (The population pH effect is tested by
    test_loose_seal_worse_than_baseline.)
    """
    kw = dict(
        storage_months=6, local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, storage_T_C=-20.0, albumin_g_dL=5.5,
        cryo_purity=1.0, glass_density=0.5, fill_volume_mL=7.5,
    )
    d_high = compute_vial_deficit_scenario(ph_storage=8.0, thaw_h=10e-6, **kw)
    d_low  = compute_vial_deficit_scenario(ph_storage=7.81, thaw_h=10e-6, **kw)
    assert abs(d_high - d_low) < 0.002

def test_vortex_thaw_reduces_deficit():
    """Smaller thaw boundary layer → faster dissolution → lower 60-min deficit."""
    kw = dict(
        storage_months=6, local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, storage_T_C=-20.0, albumin_g_dL=5.5,
        cryo_purity=1.0, glass_density=0.5, fill_volume_mL=7.5,
        ph_storage=7.81,
    )
    d_q = compute_vial_deficit_scenario(thaw_h=10e-6, **kw)
    d_v = compute_vial_deficit_scenario(thaw_h=2e-6,  **kw)
    assert d_v < d_q


# ── combined_plus and extended-mixing reference ──────────────────────────────

def test_combined_plus_at_least_as_good_as_combined_at_12mo(stats):
    """Double-pulse vortex + 90-min window is at least as effective as combined."""
    m_combined      = stats["+combined"]["mean_deficit_12mo_pct"]
    m_combined_plus = stats["+combined_plus"]["mean_deficit_12mo_pct"]
    assert m_combined_plus <= m_combined + 1e-9, (
        f"combined_plus mean {m_combined_plus} should be ≤ combined {m_combined}"
    )

def test_combined_plus_eliminates_measurable_deficit_at_12mo(stats):
    """Combined+ (mixing) brings the fraction with a measurable deficit to <5%."""
    f = stats["+combined_plus"]["frac_with_deficit_12mo"]
    assert f < 0.05, (
        f"Combined+: {f:.1%} vials with measurable deficit at 12 months; expected <5%"
    )

def test_deep_freeze_prevents_deficit(stats):
    """
    Deep-frozen / vitrified storage (≤ −80°C) arrests nucleation, so the
    precipitate never forms and even a standard quiescent thaw is clean. This is
    the root-cause PREVENTION lever (vs the mixing-step neutralization).
    """
    f6  = stats["+deep_freeze"]["frac_with_deficit_6mo"]
    f12 = stats["+deep_freeze"]["frac_with_deficit_12mo"]
    assert f6  < 0.02 and f12 < 0.02, (
        f"Deep-freeze should prevent the deficit; got {f6:.1%}/{f12:.1%}"
    )

def test_colder_storage_reduces_affected_fraction(stats):
    """Deep-freeze must leave far fewer affected vials than the −20°C baseline."""
    assert stats["+deep_freeze"]["frac_with_deficit_6mo"] < \
           stats["baseline"]["frac_with_deficit_6mo"]

def test_extended_mixing_near_zero_at_both_timepoints(stats):
    """
    Extended mixing / cold equilibration brings the fraction with a measurable
    deficit to <2% at both 6 and 12 months — consistent with the Seeker's report
    that mixing/standing reverses the deficit. Reproduced with no fitted parameter.
    """
    f6  = stats["+extended_mixing"]["frac_with_deficit_6mo"]
    f12 = stats["+extended_mixing"]["frac_with_deficit_12mo"]
    assert f6  < 0.02, f"Extended mixing: {f6:.1%} at 6mo; expected <2%"
    assert f12 < 0.02, f"Extended mixing: {f12:.1%} at 12mo; expected <2%"

def test_longer_thaw_reduces_deficit():
    """Longer thaw window (90 min vs 60 min) → more dissolution → lower deficit."""
    kw = dict(
        storage_months=12, local_k_factor=1.0, nucleation_delay_days=0.0,
        particle_size_nm=50.0, storage_T_C=-20.0, albumin_g_dL=5.5,
        cryo_purity=1.0, glass_density=0.5, fill_volume_mL=7.5,
        ph_storage=7.81, thaw_h=2e-6,
    )
    d_60  = compute_vial_deficit_scenario(thaw_min=60.0,  **kw)
    d_90  = compute_vial_deficit_scenario(thaw_min=90.0,  **kw)
    d_48h = compute_vial_deficit_scenario(thaw_min=2880.0, d_factor=0.60, **kw)
    assert d_90  < d_60,  "90-min thaw should give lower deficit than 60-min"
    assert d_48h < d_90,  "48-h cold soak should give lower deficit than 90-min"
