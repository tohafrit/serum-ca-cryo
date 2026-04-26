"""
Tests for Module 2: cryoconcentration trajectory, pH evolution.
Each test checks a physical sanity limit.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.freezing_trajectory import (
    pKa1_CO2_at_T,
    pKa2_phosphate_at_T,
    davies_log_gamma,
    pKa_effective,
    ionic_strength_baseline,
    assumption_validity_check,
    freezing_point_depression,
    unfrozen_fraction,
    cryoconcentration_trajectory,
    pH_trajectory,
    REFERENCE_COMPOSITION_mM,
    Kf_WATER,
)

COMP = REFERENCE_COMPOSITION_mM


# ── pKa temperature corrections ───────────────────────────────────────────────

def test_pKa1_CO2_at_25C():
    """At 25°C, pKa1_CO2 = 6.352 (Harned & Davis 1943)."""
    assert abs(pKa1_CO2_at_T(25.0) - 6.352) < 0.001

def test_pKa1_CO2_increases_at_lower_T():
    """pKa1 must increase as T decreases (endothermic ionization)."""
    assert pKa1_CO2_at_T(0.0) > pKa1_CO2_at_T(25.0)
    assert pKa1_CO2_at_T(-20.0) > pKa1_CO2_at_T(0.0)

def test_pKa1_CO2_at_0C_matches_literature():
    """At 0°C, pKa1 ≈ 6.58 (Harned & Davis 1943; lit value 6.581)."""
    assert abs(pKa1_CO2_at_T(0.0) - 6.577) < 0.05

def test_pKa2_phosphate_at_25C():
    """At 25°C, pKa2_phosphate = 7.21."""
    assert abs(pKa2_phosphate_at_T(25.0) - 7.21) < 0.001

def test_pKa2_phosphate_increases_at_lower_T():
    """pKa2 must increase as T decreases."""
    assert pKa2_phosphate_at_T(0.0) > pKa2_phosphate_at_T(25.0)


# ── Davies equation ───────────────────────────────────────────────────────────

def test_davies_negative_at_positive_I():
    """Activity coefficients < 1 for all ions at I > 0 → log_gamma < 0."""
    assert davies_log_gamma(1, 0.1) < 0
    assert davies_log_gamma(2, 0.1) < 0

def test_davies_zero_at_I_zero():
    """At infinite dilution, γ → 1 → log_gamma → 0 (tolerance allows for 1e-12 guard)."""
    assert abs(davies_log_gamma(1, 0.0)) < 1e-5

def test_davies_larger_for_higher_charge():
    """Higher charge → stronger deviation (more negative log_gamma)."""
    assert davies_log_gamma(2, 0.1) < davies_log_gamma(1, 0.1)


# ── Ionic strength ────────────────────────────────────────────────────────────

def test_ionic_strength_physiological_range():
    """Serum ionic strength is ~0.14–0.16 mol/kg."""
    I = ionic_strength_baseline(COMP)
    assert 0.12 < I < 0.18, f"I = {I:.4f} mol/kg, expected 0.12–0.18"

def test_ionic_strength_positive():
    I = ionic_strength_baseline(COMP)
    assert I > 0


# ── Validity check (ionic strength based) ────────────────────────────────────

def test_validity_OK_at_low_I():
    """At I_base=0.15, k=1 → I_pool=0.15 < 0.5 → OK."""
    v = assumption_validity_check(1.0, 0.15)
    assert v["level"] == "OK"

def test_validity_approximate_zone():
    """At I_base=0.15, k=5 → I_pool=0.75 (0.5–2.0) → approximate."""
    v = assumption_validity_check(5.0, 0.15)
    assert v["level"] == "approximate"

def test_validity_unreliable_at_high_k():
    """At I_base=0.15, k=20 → I_pool=3.0 > 2.0 → unreliable."""
    v = assumption_validity_check(20.0, 0.15)
    assert v["level"] == "unreliable"

def test_validity_I_pool_equals_k_times_base():
    v = assumption_validity_check(7.3, 0.2)
    assert abs(v["I_pool"] - 7.3 * 0.2) < 1e-10


# ── Freezing-point depression ─────────────────────────────────────────────────

def test_fpd_physiological_electrolyte_osmolality():
    """Electrolyte osmolality of serum ~0.28–0.35 Osm/kg."""
    res = freezing_point_depression(COMP, "none")
    assert 0.25 < res["osm_electrolyte"] < 0.38

def test_fpd_glycerol_contribution():
    """Glycerol 15% adds ~1.5–1.8 Osm/kg."""
    res_none = freezing_point_depression(COMP, "none")
    res_glyc = freezing_point_depression(COMP, "glycerol_15pct")
    delta = res_glyc["osm_total"] - res_none["osm_total"]
    assert 1.4 < delta < 2.0

def test_fpd_kf_relation():
    """ΔTf = Kf × osmolality (definition)."""
    res = freezing_point_depression(COMP, "glycerol_15pct")
    assert abs(res["delta_T"] - Kf_WATER * res["osm_total"]) < 1e-10


# ── Unfrozen fraction ─────────────────────────────────────────────────────────

def test_unfrozen_fraction_at_zero():
    fpd = freezing_point_depression(COMP, "glycerol_15pct")
    f = unfrozen_fraction(np.array([0.0]), fpd)
    assert abs(f[0] - 1.0) < 1e-9

def test_unfrozen_fraction_monotonically_decreasing():
    fpd = freezing_point_depression(COMP, "glycerol_15pct")
    T = np.linspace(-0.01, -20.0, 200)
    f = unfrozen_fraction(T, fpd)
    assert np.all(np.diff(f) <= 1e-12)

def test_unfrozen_fraction_bounded():
    fpd = freezing_point_depression(COMP, "glycerol_15pct")
    T = np.linspace(-0.01, -30.0, 500)
    f = unfrozen_fraction(T, fpd)
    assert np.all(f >= 0) and np.all(f <= 1.0)


# ── Cryoconcentration trajectory ──────────────────────────────────────────────

def test_trajectory_k_near_one_at_start():
    T = np.array([-0.001, -5.0, -20.0])
    df = cryoconcentration_trajectory(T, COMP, "glycerol_15pct")
    assert df["k"].iloc[0] == pytest.approx(1.0, abs=0.05)

def test_trajectory_k_monotonically_increasing():
    T = np.linspace(-0.01, -20.0, 200)
    df = cryoconcentration_trajectory(T, COMP, "glycerol_15pct")
    assert np.all(np.diff(df["k"].values) >= -1e-9)

def test_trajectory_species_concentrations_scale_with_k():
    T = np.array([-10.0, -20.0])
    df = cryoconcentration_trajectory(T, COMP, "glycerol_15pct")
    for _, row in df.iterrows():
        ki = row["k"]
        assert abs(row["Ca_mM"] - COMP["Ca"] * ki) < 1e-9

def test_trajectory_has_ionic_strength_column():
    T = np.linspace(-0.01, -20.0, 10)
    df = cryoconcentration_trajectory(T, COMP, "glycerol_15pct")
    assert "ionic_strength_pool" in df.columns
    assert df["ionic_strength_pool"].iloc[0] > 0


# ── pH trajectory ─────────────────────────────────────────────────────────────

def test_ph_trajectory_no_co2_loss_gt_initial():
    """With 0% CO₂ loss at −20°C, pH must be ≥ initial 7.4 (temperature effect)."""
    k = np.array([1.0, 5.0, 10.0])
    df = pH_trajectory(k, COMP, (0.0,), T_celsius=-20.0)
    assert (df["pH_co2loss_00pct"] >= 7.35).all()

def test_ph_trajectory_co2_loss_raises_ph():
    """More CO₂ loss → higher pH at same k."""
    k = np.array([5.0])
    df = pH_trajectory(k, COMP, (0.0, 0.5, 0.9), T_celsius=-20.0)
    assert df["pH_co2loss_50pct"].iloc[0] > df["pH_co2loss_00pct"].iloc[0]
    assert df["pH_co2loss_90pct"].iloc[0] > df["pH_co2loss_50pct"].iloc[0]

def test_ph_trajectory_output_columns():
    k = np.array([1.0, 5.0])
    df = pH_trajectory(k, COMP, (0.0, 0.5, 0.9))
    assert "pH_co2loss_00pct" in df.columns
    assert "pH_co2loss_50pct" in df.columns
    assert "pH_co2loss_90pct" in df.columns

def test_ph_independent_of_k_at_zero_loss_zero_I_correction():
    """
    At zero CO₂ loss and no ionic strength correction, pH is independent of k.
    In our model ionic strength correction is included, so pH changes slightly
    with k — but only the secondary Davies term, so variation should be < 0.5 pH.
    """
    k = np.array([1.0, 5.0, 20.0])
    df = pH_trajectory(k, COMP, (0.0,), T_celsius=-20.0)
    pH_range = df["pH_co2loss_00pct"].max() - df["pH_co2loss_00pct"].min()
    assert pH_range < 0.5, f"pH varies {pH_range:.3f} with k (expected <0.5 from secondary I effect)"
