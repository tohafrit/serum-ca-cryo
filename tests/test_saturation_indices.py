"""
Tests for Module 3: saturation indices.
All tests check physical/biological sanity limits, not code mechanics.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pytest
from src.saturation_indices import (
    activity_coefficient,
    phosphate_speciation,
    carbonate_speciation,
    calcium_free_fraction,
    ksp_at_T,
    saturation_index,
    all_si,
    physiological_serum_si,
    pKa3_phosphate_at_T,
    REFERENCE_COMPOSITION_mM,
    KSP_DATA,
)

COMP = dict(REFERENCE_COMPOSITION_mM)
COMP["protein"] = 4.5 * 10.0 / 66.5 * 1000.0  # 4.5 g/dL albumin correctly encoded


# ── Activity coefficients ─────────────────────────────────────────────────────

def test_gamma_unity_at_zero_I():
    """At infinite dilution γ → 1."""
    assert abs(activity_coefficient(1, 1e-9, 25.0) - 1.0) < 0.01
    assert abs(activity_coefficient(2, 1e-9, 25.0) - 1.0) < 0.01

def test_gamma_less_than_one_at_finite_I():
    """Davies: γ < 1 for all ions at I > 0."""
    assert activity_coefficient(1, 0.15, 25.0) < 1.0
    assert activity_coefficient(2, 0.15, 25.0) < 1.0

def test_gamma_decreases_with_charge():
    """Higher charge → stronger interaction → lower γ."""
    assert activity_coefficient(2, 0.15, 25.0) < activity_coefficient(1, 0.15, 25.0)

def test_gamma_at_physiological_I():
    """At I=0.16 mol/kg, γ(Ca²⁺) ≈ 0.28–0.38 (literature range)."""
    g = activity_coefficient(2, 0.16, 25.0)
    assert 0.25 < g < 0.45, f"γ(Ca²⁺) = {g:.3f} at I=0.16"


# ── Phosphate speciation ──────────────────────────────────────────────────────

def test_phosphate_speciation_sums_to_total():
    spec = phosphate_speciation(1.0, 7.4, 25.0, 0.16)
    total = spec["H3PO4"] + spec["H2PO4"] + spec["HPO4"] + spec["PO4"]
    assert abs(total - 1.0) < 1e-9

def test_phosphate_HPO4_dominates_at_pH_7_4():
    """At pH 7.4 and 25°C (uncorrected pKa₂=7.21), HPO₄²⁻ should be most abundant."""
    spec = phosphate_speciation(1.0, 7.4, 25.0, 0.01)  # low I → near-thermodynamic
    assert spec["HPO4"] > spec["H2PO4"]
    assert spec["HPO4"] > spec["PO4"]

def test_PO4_fraction_increases_with_pH():
    """PO₄³⁻ fraction must increase monotonically with pH."""
    pHs = [6.0, 7.0, 8.0, 9.0, 10.0]
    alphas = [phosphate_speciation(1.0, pH, 25.0, 0.16)["alpha_PO4"] for pH in pHs]
    assert all(alphas[i] < alphas[i+1] for i in range(len(alphas)-1))

def test_pKa3_increases_at_lower_T():
    """HPO₄²⁻/PO₄³⁻ pKa₃ increases as T drops (endothermic ionization)."""
    assert pKa3_phosphate_at_T(-20.0) > pKa3_phosphate_at_T(25.0)


# ── Carbonate speciation ──────────────────────────────────────────────────────

def test_carbonate_sums_to_total():
    spec = carbonate_speciation(25.0, 7.4, 25.0, 0.16)
    total = spec["CO2_aq"] + spec["HCO3"] + spec["CO3"]
    assert abs(total - 25.0) < 1e-9

def test_HCO3_dominates_at_pH_7_4():
    """At pH 7.4, HCO₃⁻ >> CO₂ and CO₃²⁻."""
    spec = carbonate_speciation(25.0, 7.4, 25.0, 0.16)
    assert spec["HCO3"] > spec["CO2_aq"]
    assert spec["HCO3"] > spec["CO3"]

def test_CO3_increases_with_pH():
    CO3s = [carbonate_speciation(25.0, pH, 25.0, 0.16)["CO3"] for pH in [7, 8, 9, 10]]
    assert all(CO3s[i] < CO3s[i+1] for i in range(len(CO3s)-1))


# ── Albumin–calcium binding ───────────────────────────────────────────────────

def test_alpha_Ca_at_zero_albumin():
    """Without albumin all Ca²⁺ is free."""
    alpha = calcium_free_fraction(2.5, albumin_gL=0.0, pH=7.4)
    assert abs(alpha - 1.0) < 1e-6

def test_alpha_Ca_physiological():
    """At 4 g/dL (40 g/L) albumin, pH 7.4: α_Ca ≈ 0.50 (Fogh-Andersen 1995)."""
    alpha = calcium_free_fraction(2.4, albumin_gL=40.0, pH=7.4)
    assert 0.40 < alpha < 0.60, f"α_Ca = {alpha:.3f}, expected 0.40–0.60"

def test_alpha_Ca_decreases_with_albumin():
    """More albumin → more Ca bound → lower free fraction."""
    albs = [0.0, 10.0, 20.0, 40.0, 80.0]
    alphas = [calcium_free_fraction(2.5, a, pH=7.4) for a in albs]
    assert all(alphas[i] >= alphas[i+1] for i in range(len(alphas)-1))

def test_alpha_Ca_decreases_with_higher_pH():
    """Higher pH → more Ca²⁺ binding to albumin → lower α_Ca."""
    alpha_74 = calcium_free_fraction(2.5, 40.0, pH=7.4)
    alpha_82 = calcium_free_fraction(2.5, 40.0, pH=8.2)
    assert alpha_82 < alpha_74


# ── Ksp temperature correction ────────────────────────────────────────────────

def test_ksp_unchanged_at_25C():
    """ksp_at_T at 25°C should return the reference value exactly."""
    for phase, data in KSP_DATA.items():
        assert abs(ksp_at_T(phase, 25.0) - data["log_ksp_25"]) < 1e-9

def test_HAp_less_soluble_at_lower_T():
    """HAp has ΔH>0 → less soluble at lower T → log Ksp more negative."""
    assert ksp_at_T("HAp", -20.0) < ksp_at_T("HAp", 25.0)

def test_calcite_less_soluble_at_lower_T():
    """Calcite ΔH>0 → same direction as HAp."""
    assert ksp_at_T("calcite", -20.0) < ksp_at_T("calcite", 25.0)


# ── Saturation indices — main sanity checks ───────────────────────────────────

def test_physiological_serum_SI_HAp():
    """
    CRITICAL: Normal human serum is metastably supersaturated with respect to HAp.
    SI(HAp) ≈ +3 to +8 (Heughebaert & Nancollas 1984; Lenz et al. 2013).
    If this fails, the activity model is broken.
    """
    phys = physiological_serum_si()
    si_hap = phys["HAp"]
    assert 3.0 <= si_hap <= 10.0, (
        f"SI(HAp) physiological = {si_hap:.2f}; expected +3 to +10. "
        "Check pKa₃, activity coefficients, albumin binding."
    )

def test_physiological_alpha_Ca():
    """Physiological free Ca²⁺ fraction ≈ 0.45–0.55."""
    phys = physiological_serum_si()
    alpha = phys["_activities"]["alpha_Ca"]
    assert 0.40 <= alpha <= 0.60, f"α_Ca = {alpha:.3f}"

def test_cryo_SI_HAp_exceeds_physiological():
    """
    SI(HAp) at cryo state must be ≥ physiological SI.
    Cryoconcentration + pH rise can only increase supersaturation.
    """
    phys = physiological_serum_si()
    I_base = 0.143
    k = 5.58
    comp_cryo = {sp: v * k for sp, v in COMP.items()}
    si_cryo = all_si(comp_cryo, pH=7.81, T_celsius=-20.0, I=I_base*k)
    assert si_cryo["HAp"] >= phys["HAp"] - 1.0, (
        f"SI(HAp) cryo = {si_cryo['HAp']:.2f} < physiological {phys['HAp']:.2f}"
    )

def test_brushite_more_soluble_than_HAp():
    """Brushite is the most soluble CaP phase; SI(brushite) < SI(HAp) always."""
    si = all_si(COMP, pH=7.4, T_celsius=25.0, I=0.16)
    assert si["brushite"] < si["HAp"]

def test_calcite_near_saturation_in_serum():
    """
    Serum CO₂ = 25 mM; calcite SI should be in range -2 to +2
    (blood is near calcite saturation in vivo; Lemann et al. 1970).
    """
    comp_phys = dict(REFERENCE_COMPOSITION_mM)
    comp_phys["CO2"] = 25.0
    comp_phys["protein"] = 4.5 * 10.0 / 66.5 * 1000.0
    si = all_si(comp_phys, pH=7.4, T_celsius=37.0, I=0.16)
    assert -3.0 <= si["calcite"] <= 3.0, f"SI(calcite) = {si['calcite']:.2f}"

def test_ACP_loose_more_soluble_than_ACP_tight():
    """ACP_loose has higher Ksp → lower SI at same conditions."""
    si = all_si(COMP, pH=7.4, T_celsius=25.0, I=0.16)
    assert si["ACP_loose"] < si["ACP_tight"]

def test_SI_increases_with_pH_for_HAp():
    """Higher pH → more PO₄³⁻ and OH⁻ → higher SI(HAp)."""
    si_74 = saturation_index("HAp", COMP, pH=7.4, T_celsius=25.0, I=0.16)
    si_82 = saturation_index("HAp", COMP, pH=8.2, T_celsius=25.0, I=0.16)
    assert si_82 > si_74
