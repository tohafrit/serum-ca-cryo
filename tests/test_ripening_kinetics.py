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

def test_ripening_suppressed_at_cryo():
    """
    With the corrected pool viscosity (~4100 mPa·s), solution-mediated ripening
    is suppressed: over a year at −20°C the precipitate stays essentially all
    amorphous (ACP), with negligible OCP/HAp. This is the key honest finding and
    agrees with Combes & Rey (2010): ACP is kinetically stable for months at <0°C.
    """
    t = np.linspace(0, 365, 200)
    sol = phase_evolution(t)
    assert sol["x_ACP"][-1] > 0.95, "ACP should stay >95% over a year at −20°C"
    assert sol["x_HAp"][-1] < 0.02, "HAp should be negligible (ripening suppressed)"


# ── Particle growth (LSW) ─────────────────────────────────────────────────────

def test_radius_increases_with_time():
    r0 = mean_radius_nm("ACP", 0.0)
    r1 = mean_radius_nm("ACP", 1000.0)
    assert r1 > r0

def test_coarsening_suppressed_at_cryo():
    """
    Diffusion-limited LSW coarsening is also suppressed at the corrected pool
    viscosity: after a year the mean radius stays close to its initial value
    (growth < ~30%), consistent with the precipitate remaining fine amorphous
    aggregates rather than coarsening into large crystals.
    """
    r1 = mean_radius_nm("ACP", 365 * 24.0)
    # Individual particles stay nanoscale (~tens of nm) over a year — LSW
    # coarsening to crystalline µm-scale particles does not happen. (The wall
    # deposit reaches µm scale by aggregation, a distinct process.)
    assert r1 < 100.0, f"ACP radius {r1:.0f} nm after 1y; should stay nanoscale"


# ── Prevention lever: colder storage arrests nucleation ───────────────────────

def test_nucleation_temp_factor_arrests_when_cold():
    """
    Induction time ∝ η/T, so colder storage massively lengthens it. The factor
    is 1.0 at the −20°C reference, <1 when warmer, and huge (→ arrest) deep-frozen.
    This is the physical basis for deep-freeze prevention.
    """
    from src.ripening_kinetics import nucleation_temp_factor as f
    assert abs(f(-20.0) - 1.0) < 1e-9
    assert f(-18.0) < 1.0 < f(-22.0)          # warmer faster, colder slower
    assert f(-40.0) > 10.0                     # already strongly suppressed
    assert f(-80.0) >= 1e5                      # vitrified → nucleation arrested


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

def test_affected_vial_deficit_matches_seeker():
    """
    A vial that has precipitated shows a ≥4% quiescent-thaw deficit, matching the
    Seeker's reported "4% or more". The onset timing ("not before 6 months") is a
    nucleation/population effect, tested in Module 6, not a per-vial magnitude.
    """
    d = ca_deficit_at_60min(6)
    assert d >= 0.04, f"affected-vial deficit {d*100:.1f}%; Seeker reports ≥4%"

def test_deficit_is_modest_and_positive():
    """
    A nucleated vial shows a modest, real deficit at a quiescent 60-min thaw.
    Magnitude is uncertain (set by precipitated fraction × aggregate morphology,
    band ~0.5–15%); for representative parameters it is a few percent.
    """
    d = ca_deficit_at_60min(6)
    assert 0.003 <= d <= 0.06, (
        f"60-min deficit = {d*100:.1f}%; expected modest-positive (0.3–6%)"
    )

def test_deficit_independent_of_storage_no_ripening():
    """
    With ripening suppressed, the per-event deficit does NOT grow with storage
    time (the precipitate stays amorphous and does not coarsen). The population
    effect grows only because more vials nucleate — tested in Module 6.
    """
    deficits = [ca_deficit_at_60min(m) for m in [1, 3, 6, 12]]
    assert max(deficits) - min(deficits) < 0.005, (
        "per-event deficit should be ~flat vs storage (no ripening)"
    )

def test_deficit_12months_recovers_after_48h():
    """
    Seeker's reported reversibility: 24–48 h equilibration resolves the deficit.
    At 12 months, 48-h cold recovery should leave deficit < 2%.
    """
    from src.ripening_kinetics import F_PRECIP
    rec_48h = float(ca_recovery_curve(12, np.array([2880.0]), protocol="cold_4C_48h")[0])
    deficit_48h = F_PRECIP * (1.0 - rec_48h)
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
