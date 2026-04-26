"""
Tests for Module 5: Ostwald ripening kinetics and post-thaw redissolution.
All tests check physical/biological sanity — no code mechanics.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.ripening_kinetics import (
    rate_constants,
    phase_evolution,
    ca_recovery_curve,
    ca_deficit_at_60min,
    mean_radius_nm,
    c_sat_ostwald,
)


# ── Rate constants ────────────────────────────────────────────────────────────

def test_rate_constants_slower_at_minus20_than_25C():
    """Arrhenius: all rates decrease at −20°C vs 25°C."""
    k25  = rate_constants(25.0,  7.0)
    km20 = rate_constants(-20.0, 7.0)
    for key in k25:
        assert km20[key] < k25[key], f"{key} should be slower at −20°C"

def test_rate_constants_faster_at_higher_pH():
    """ACP→OCP rate increases with pH (Christoffersen 1989)."""
    k_low  = rate_constants(-20.0, 7.0)
    k_high = rate_constants(-20.0, 8.0)
    assert k_high["k_ACP_OCP"] > k_low["k_ACP_OCP"]

def test_rate_constants_positive():
    kk = rate_constants(-20.0, 7.81)
    for key, val in kk.items():
        assert val > 0, f"{key} must be positive"


# ── Phase evolution ODE ───────────────────────────────────────────────────────

def test_phase_fractions_sum_to_one():
    t = np.linspace(0, 730, 200)
    sol = phase_evolution(t)
    total = sol["x_ACP"] + sol["x_OCP"] + sol["x_HAp"]
    assert np.allclose(total, 1.0, atol=1e-6), "Phase fractions must sum to 1"

def test_ACP_monotonically_decreasing():
    t = np.linspace(0, 730, 200)
    sol = phase_evolution(t)
    assert np.all(np.diff(sol["x_ACP"]) <= 1e-9)

def test_HAp_monotonically_increasing():
    t = np.linspace(0, 730, 200)
    sol = phase_evolution(t)
    assert np.all(np.diff(sol["x_HAp"]) >= -1e-9)

def test_initial_conditions_pure_ACP():
    sol = phase_evolution(np.array([0.0, 1.0]))
    assert abs(sol["x_ACP"][0] - 1.0) < 1e-9
    assert abs(sol["x_OCP"][0]) < 1e-9
    assert abs(sol["x_HAp"][0]) < 1e-9

def test_ACP_majority_at_1_month():
    """ACP should still dominate at 1 month (< 30 days)."""
    sol = phase_evolution(np.array([0.0, 30.0]))
    assert sol["x_ACP"][-1] > 0.7, "ACP should be >70% at 1 month"

def test_OCP_rises_then_eventually_declines():
    """OCP should peak and then (very slowly) decline as HAp takes over."""
    t = np.linspace(0, 5000, 500)
    sol = phase_evolution(t)
    peak_idx = np.argmax(sol["x_OCP"])
    # OCP peak should occur after t=0 but before the end
    assert 0 < peak_idx < len(t) - 1


# ── Particle growth (LSW) ─────────────────────────────────────────────────────

def test_radius_increases_with_time():
    r0 = mean_radius_nm("ACP", 0.0)
    r1 = mean_radius_nm("ACP", 1000.0)
    assert r1 > r0

def test_radius_ACP_larger_than_HAp_at_same_time():
    """ACP grows faster (lower surface energy, higher LSW K)."""
    assert mean_radius_nm("ACP", 1000) > mean_radius_nm("HAp", 1000)


# ── Ostwald-Freundlich solubility ─────────────────────────────────────────────

def test_c_sat_decreases_with_radius():
    """Larger particles → lower curvature → lower effective solubility."""
    c_small = c_sat_ostwald("HAp", 10.0)
    c_large = c_sat_ostwald("HAp", 500.0)
    assert c_small > c_large

def test_c_sat_HAp_below_ACP():
    """HAp is less soluble than ACP at the same particle size."""
    r = 100.0
    assert c_sat_ostwald("HAp", r) < c_sat_ostwald("ACP", r)


# ── Ca recovery curves ────────────────────────────────────────────────────────

def test_recovery_bounded_zero_to_one():
    t = np.linspace(0, 60, 100)
    rec = ca_recovery_curve(6, t, protocol="quiescent_22C")
    assert np.all(rec >= 0) and np.all(rec <= 1.0)

def test_recovery_monotonically_increasing():
    t = np.linspace(0, 60, 100)
    rec = ca_recovery_curve(6, t, protocol="quiescent_22C")
    assert np.all(np.diff(rec) >= -1e-9)

def test_recovery_at_zero_time_is_zero():
    rec = ca_recovery_curve(6, np.array([0.0]))
    assert abs(float(rec[0])) < 1e-6

def test_vortex_recovers_more_than_quiescent():
    """Smaller boundary layer → faster dissolution → higher recovery at 60 min."""
    t = np.array([60.0])
    rec_q = float(ca_recovery_curve(12, t, protocol="quiescent_22C")[0])
    rec_v = float(ca_recovery_curve(12, t, protocol="vortex_22C")[0])
    assert rec_v > rec_q, "Vortex should improve 60-min recovery"

def test_cold_protocol_approaches_full_recovery_at_48h():
    """At 2–8°C over 48 h, even 12-month samples should recover >90%."""
    rec = ca_recovery_curve(12, np.array([2880.0]), protocol="cold_4C_48h")
    assert float(rec[0]) > 0.90, "48-h cold equilibration should give >90% recovery"


# ── Ca deficit — the critical falsification tests ─────────────────────────────

def test_deficit_1month_below_threshold():
    """At 1 month storage: deficit < 2% (no observable effect)."""
    d = ca_deficit_at_60min(1)
    assert d < 0.02, f"1-month deficit = {d*100:.1f}%; expected < 2%"

def test_deficit_6months_brackets_seeker_threshold():
    """
    CRITICAL: At 6 months, 60-min deficit must be 3–8%.
    Seeker reports '≥4%' as the observed phenomenon.
    If this fails, the kinetic parameters need revision.
    """
    d = ca_deficit_at_60min(6)
    assert 0.03 <= d <= 0.08, (
        f"6-month deficit = {d*100:.1f}%; expected 3–8% "
        "(Seeker reports ≥4%). Check Arrhenius Ea, rate constants, f_precip."
    )

def test_deficit_increases_with_storage_time():
    """More months → more HAp → larger 60-min deficit."""
    deficits = [ca_deficit_at_60min(m) for m in [1, 3, 6, 12]]
    assert all(deficits[i] < deficits[i+1] for i in range(len(deficits)-1))

def test_deficit_12months_recovers_after_48h():
    """
    Seeker's reported reversibility: 24–48 h equilibration resolves the deficit.
    At 12 months, 48-h cold recovery should leave deficit < 2%.
    """
    rec_48h = float(ca_recovery_curve(12, np.array([2880.0]), protocol="cold_4C_48h")[0])
    deficit_48h = 0.88 * (1.0 - rec_48h)
    assert deficit_48h < 0.02, (
        f"48-h residual deficit = {deficit_48h*100:.1f}%; expected < 2%"
    )

def test_vortex_reduces_6month_deficit_meaningfully():
    """30-s vortex should reduce the 6-month, 60-min deficit by ≥20%."""
    d_q = ca_deficit_at_60min(6, protocol="quiescent_22C")
    d_v = ca_deficit_at_60min(6, protocol="vortex_22C")
    reduction = (d_q - d_v) / d_q
    assert reduction >= 0.20, (
        f"Vortex reduces deficit by {reduction*100:.0f}%; expected ≥20%"
    )
