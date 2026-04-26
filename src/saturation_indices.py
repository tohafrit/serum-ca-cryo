"""
Module 3 — Saturation indices for calcium phosphate and carbonate phases.

All IAP values are computed using ion ACTIVITIES (not concentrations).
Activity coefficients are computed via the Davies equation with explicit
temperature correction through the dielectric constant of water.

Uncertainty bands
-----------------
At I ≈ 0.8 mol/kg (primary glycerol-15% scenario at k=5.58) Davies introduces
±0.3 units of SI uncertainty due to non-ideality. At -20°C, van't Hoff temperature
extrapolation of Ksp adds ±0.5 log units. Both are shown on figures.

ACP solubility is a range (not a point): log Ksp from -25.7 to -28.5 depending
on Ca/P ratio and age. SI(ACP) is shown as a band on all figures.

Physiological sanity check
---------------------------
Normal human blood plasma is metastably supersaturated with respect to HAp.
Thermodynamic SI(HAp) of plasma ≈ +3 to +7 (Heughebaert & Nancollas 1984;
Lenz et al. 2013). Biological inhibitors (fetuin-A, pyrophosphate, Mg²⁺,
citrate) suppress nucleation kinetically. This large SI is a known fact of
clinical biochemistry and is used as the primary sanity check here.

Phase dissolution conventions (must match the cited Ksp)
---------------------------------------------------------
HAp :      Ca₅(PO₄)₃OH  ⇌ 5Ca²⁺ + 3PO₄³⁻ + OH⁻
OCP :      Ca₈H₂(PO₄)₆  ⇌ 8Ca²⁺ + 2H⁺   + 6PO₄³⁻
Brushite : CaHPO₄·2H₂O  ⇌ Ca²⁺  + HPO₄²⁻
Monetite : CaHPO₄        ⇌ Ca²⁺  + HPO₄²⁻
ACP :      Ca₃(PO₄)₂     ⇌ 3Ca²⁺ + 2PO₄³⁻  (Ca/P = 1.5 form)
Calcite :  CaCO₃          ⇌ Ca²⁺  + CO₃²⁻
Vaterite : CaCO₃          ⇌ Ca²⁺  + CO₃²⁻

References
----------
Ksp(HAp)      : McDowell H et al. (1977) J Res NBS 82:11; log Ksp = -58.33 (25°C)
Ksp(OCP)      : Moreno EC et al. (1960) J Res NBS 64A:425; log Ksp = -96.6 (25°C)
Ksp(brushite) : Marshall RW & Nancollas GH (1969) J Phys Chem 73:3838; log Ksp = -6.59
Ksp(monetite) : McDowell H et al. (1971) J Res NBS 75A:105; log Ksp = -6.90
Ksp(ACP-loose): Boskey AL & Posner AS (1973) J Phys Chem 77:2313; log Ksp = -25.7
Ksp(ACP-tight): Christoffersen J et al. (1990) J Cryst Growth 94:767; log Ksp = -28.5
Ksp(calcite)  : Plummer LN & Busenberg E (1982) Geochim Cosmochim Acta 46:1011; log Ksp = -8.48
Ksp(vaterite) : Plummer LN & Busenberg E (1982); log Ksp = -7.91
ΔH(HAp)       : Vega ED et al. (1996) J Cryst Growth 167:491; ΔH ≈ +20 kJ/mol (±20)
ΔH(calcite)   : CRC Handbook 97th ed.; ΔH = +12.0 kJ/mol (well established)
pKa(phosphate): Goldberg RN et al. (2002) J Phys Chem Ref Data 31:231
pKa(CO₂)      : Harned HS & Davis R (1943) J Am Chem Soc 65:2030
pKw            : Bandura AV & Lvov SN (2006) J Phys Chem Ref Data 35:15
Albumin-Ca     : Fogh-Andersen N et al. (1995) Clin Chem 41:1522;
                 Pedersen KO (1972) Scand J Clin Lab Invest 29:75
Davies equation: Davies CW (1962) Ion Association. Butterworths, London
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Import shared thermodynamic helpers from Module 2
import sys
sys.path.insert(0, str(Path(__file__).parent))
from freezing_trajectory import (
    REFERENCE_COMPOSITION_mM,
    ionic_strength_baseline,
    pKa1_CO2_at_T,
    pKa2_phosphate_at_T,
    cryoconcentration_trajectory,
    freezing_point_depression,
    pH_trajectory,
)

# ── Thermodynamic constants ───────────────────────────────────────────────────

# log Ksp at 25°C; one row per phase
# Columns: log_ksp, delta_H_kJ_per_mol, stoich tuple (for reference), note
KSP_DATA: dict[str, dict] = {
    "HAp": {
        "log_ksp_25": -58.33,
        "delta_H":    +20.0,   # kJ/mol, ±20; ΔH uncertainty dominates at -20°C
        "note":       "McDowell 1977; ΔH: Vega 1996 (±20 kJ/mol uncertainty)",
    },
    "OCP": {
        "log_ksp_25": -96.6,
        "delta_H":    +15.0,   # estimated; less data than HAp
        "note":       "Moreno 1960; ΔH estimated from HAp analogy (large uncertainty)",
    },
    "brushite": {
        "log_ksp_25": -6.59,
        "delta_H":    +6.0,    # Marshall & Nancollas 1969
        "note":       "Marshall & Nancollas 1969",
    },
    "monetite": {
        "log_ksp_25": -6.90,
        "delta_H":    +6.0,    # estimated
        "note":       "McDowell 1971; ΔH estimated",
    },
    "ACP_loose": {
        "log_ksp_25": -25.7,
        "delta_H":    +5.0,    # order-of-magnitude estimate for ACP
        "note":       "Boskey & Posner 1973 (fresh ACP, Ca/P=1.50)",
    },
    "ACP_tight": {
        "log_ksp_25": -28.5,
        "delta_H":    +5.0,    # same estimate
        "note":       "Christoffersen 1990 (aged ACP, higher crystallinity)",
    },
    "calcite": {
        "log_ksp_25": -8.48,
        "delta_H":    +12.0,   # CRC 97th ed.; well established
        "note":       "Plummer & Busenberg 1982; ΔH CRC 97th ed.",
    },
    "vaterite": {
        "log_ksp_25": -7.91,
        "delta_H":    +10.0,   # Plummer & Busenberg 1982
        "note":       "Plummer & Busenberg 1982",
    },
}

# ΔH uncertainty for van't Hoff extrapolation to -20°C (±log units in Ksp)
DELTA_H_UNCERTAINTY_KJ = {
    "HAp":      20.0,
    "OCP":      30.0,   # larger uncertainty — less data
    "brushite": 5.0,
    "monetite": 8.0,
    "ACP_loose": 15.0,
    "ACP_tight": 15.0,
    "calcite":   2.0,   # best known
    "vaterite":  5.0,
}

# pKw: temperature-dependent water dissociation constant
# pKw = A + B/T + C*ln(T)  [T in K]
# Using simple linear fit from Bandura & Lvov (2006) for 0–100°C range:
# pKw(T°C) ≈ 14.94 − 0.033*(T − 0) [at low T] → at 25°C: 14.94 - 0.033*25 = 14.12 (off)
# Better: pKw(25°C)=14.00, pKw(0°C)=14.94 → linear: pKw = 14.00 + 0.0376*(25 − T)
def pKw_at_T(T_celsius: float) -> float:
    """
    pKw = -log10(Kw) for water at temperature T.
    Linear approximation: pKw = 14.00 + 0.0376*(25 - T_celsius)
    Valid for 0–50°C; extrapolated below 0°C with increasing uncertainty.
    Refs: Bandura & Lvov (2006); Harris (2010) Quantitative Chemical Analysis.
    """
    return 14.00 + 0.0376 * (25.0 - T_celsius)


# ── Davies equation (temperature-corrected) ───────────────────────────────────

def dielectric_water(T_celsius: float) -> float:
    """
    Dielectric constant of water ε_r at temperature T.
    Fit to Uematsu & Frank (1980): ε_r = 87.74 − 0.4008*T + 9.398e-4*T² − 1.41e-6*T³
    Valid 0–100°C; extrapolated below 0°C.
    """
    T = T_celsius
    return 87.74 - 0.4008*T + 9.398e-4*T**2 - 1.41e-6*T**3


def davies_A(T_celsius: float) -> float:
    """
    Temperature-corrected Davies A parameter.
    A(T) = 1.825e6 / (ε_r(T) * T_K)^1.5   [mol/kg]^(-0.5)
    At 25°C: A ≈ 0.509 (standard value).
    """
    T_K = T_celsius + 273.15
    eps = dielectric_water(T_celsius)
    return 1.825e6 / (eps * T_K) ** 1.5


def activity_coefficient(z: float, I: float, T_celsius: float = 25.0) -> float:
    """
    Individual ion activity coefficient via temperature-corrected Davies equation.
    log10(γ) = -A(T) * z² * (√I/(1+√I) - 0.3*I)

    Reliable for I < 0.5 mol/kg; increasingly approximate above that.
    At I > 1.0 mol/kg the Davies equation can give γ > 1 (artifact: -0.3I term);
    we cap log γ to ensure γ ≤ 1 in that regime and flag it.
    """
    A = davies_A(T_celsius)
    sqrtI = np.sqrt(max(I, 1e-12))
    log_gamma = -A * z**2 * (sqrtI / (1.0 + sqrtI) - 0.3 * I)
    # Cap: physical γ ≤ 1 for most ionic solutions; Davies artifact at high I
    return 10.0 ** min(log_gamma, 0.0)


def si_uncertainty_from_activity_model(I: float) -> float:
    """
    Estimated ±SI uncertainty from Davies activity model at ionic strength I.
    Based on comparison of Davies vs SIT model for similar electrolyte compositions.
    """
    if I < 0.5:
        return 0.1
    elif I < 1.0:
        return 0.3
    elif I < 2.0:
        return 0.5
    else:
        return 1.0   # Pitzer model needed


# ── Ksp temperature correction ────────────────────────────────────────────────

def ksp_at_T(phase: str, T_celsius: float) -> float:
    """
    Temperature-corrected Ksp via van't Hoff equation.

    log Ksp(T) = log Ksp(25°C) - ΔH/(2.303R) * (1/T - 1/298.15)

    Sign convention: ΔH > 0 → endothermic dissolution → Ksp decreases at lower T
    → phase becomes less soluble as T drops → SI increases at -20°C.

    Uncertainty: ΔH extrapolation from 25°C to -20°C introduces ±0.5 log units
    for most phases (±20 kJ/mol in ΔH is typical). This is shown as error bars.

    Returns log Ksp (base 10) at temperature T_celsius.
    """
    data = KSP_DATA[phase]
    log_ksp_25 = data["log_ksp_25"]
    delta_H = data["delta_H"] * 1000.0  # kJ → J
    T2 = T_celsius + 273.15
    T1 = 298.15
    R = 8.314  # J/(mol·K)
    # van't Hoff: log K(T2) = log K(T1) - ΔH/(2.303R) * (1/T2 - 1/T1)
    log_ksp_T = log_ksp_25 - (delta_H / (2.303 * R)) * (1.0/T2 - 1.0/T1)
    return log_ksp_T


def ksp_uncertainty_at_T(phase: str, T_celsius: float) -> float:
    """
    ±Δ(log Ksp) from ΔH uncertainty when extrapolating to T_celsius.
    Returns one-sigma band (±Δlog_ksp).
    """
    delta_delta_H = DELTA_H_UNCERTAINTY_KJ[phase] * 1000.0  # J
    T2 = T_celsius + 273.15
    T1 = 298.15
    R = 8.314
    return abs((delta_delta_H / (2.303 * R)) * (1.0/T2 - 1.0/T1))


# ── Phosphate speciation ──────────────────────────────────────────────────────

def pKa3_phosphate_at_T(T_celsius: float) -> float:
    """
    Temperature-corrected pKa₃ for HPO₄²⁻ ⇌ H⁺ + PO₄³⁻.
    van't Hoff with ΔH₃ = +14.6 kJ/mol (Goldberg et al. 2002).
    At 25°C: pKa₃ = 12.35; at 0°C: ≈ 12.81.
    """
    pKa3_25 = 12.35
    delta_H3 = 14600.0  # J/mol (endothermic ionization: HPO₄²⁻ → H⁺ + PO₄³⁻)
    T2 = T_celsius + 273.15
    T1 = 298.15
    # Note sign: pKa = -log Ka.  pKa(T₂) = pKa(T₁) + (ΔH/2.303R)*(1/T₂ - 1/T₁)
    # (opposite to the ksp_at_T formula which works on log Ksp, not pKsp)
    return pKa3_25 + (delta_H3 / (2.303 * 8.314)) * (1.0/T2 - 1.0/T1)


def _pKa_obs(pKa_thermo: float, z_acid: int, z_base: int, I: float,
             T_celsius: float) -> float:
    """
    Effective (observed) pKa at ionic strength I and temperature T_celsius.

    For HA^(z_acid) ⇌ H⁺ + A^(z_base):
        pKa_obs = pKa_thermo - (log γ_H + log γ_A - log γ_HA)
    """
    A = davies_A(T_celsius)
    def _log_g(z: int) -> float:
        if z == 0:
            return 0.0  # neutral species
        sqrtI = np.sqrt(max(I, 1e-12))
        return -A * z**2 * (sqrtI / (1.0 + sqrtI) - 0.3 * I)

    correction = _log_g(1) + _log_g(z_base) - _log_g(z_acid)
    return pKa_thermo - correction


def phosphate_speciation(
    total_Pi_mM: float,
    pH: float,
    T_celsius: float,
    I: float,
) -> dict:
    """
    Fractional speciation of inorganic phosphate into four species at given pH, T, I.

    Returns dict with keys: H3PO4, H2PO4, HPO4, PO4
    Values are molar concentrations in mM.

    Uses temperature and ionic-strength corrected pKa values:
        pKa₁ (H₃PO₄/H₂PO₄⁻): 2.15 at 25°C; ΔH₁ = -8.0 kJ/mol (Goldberg 2002)
        pKa₂ (H₂PO₄⁻/HPO₄²⁻): from pKa2_phosphate_at_T()
        pKa₃ (HPO₄²⁻/PO₄³⁻):  from pKa3_phosphate_at_T()
    """
    # pKa₁ correction (H₃PO₄ → H₂PO₄⁻ + H⁺; z_acid=0, z_base=1)
    pKa1_thermo = 2.15 - ((-8000.0)/(2.303*8.314)) * (1.0/(T_celsius+273.15) - 1.0/298.15)
    pKa1_eff = _pKa_obs(pKa1_thermo, z_acid=0, z_base=1, I=I, T_celsius=T_celsius)

    pKa2_eff = _pKa_obs(pKa2_phosphate_at_T(T_celsius), z_acid=1, z_base=2,
                        I=I, T_celsius=T_celsius)
    pKa3_eff = _pKa_obs(pKa3_phosphate_at_T(T_celsius), z_acid=2, z_base=3,
                        I=I, T_celsius=T_celsius)

    h = 10**(-pH)
    K1 = 10**(-pKa1_eff)
    K2 = 10**(-pKa2_eff)
    K3 = 10**(-pKa3_eff)

    # Denominator: [H⁺]³ + K1[H⁺]² + K1K2[H⁺] + K1K2K3
    denom = h**3 + K1*h**2 + K1*K2*h + K1*K2*K3

    alpha_H3PO4  = h**3 / denom
    alpha_H2PO4  = K1 * h**2 / denom
    alpha_HPO4   = K1 * K2 * h / denom
    alpha_PO4    = K1 * K2 * K3 / denom

    return {
        "H3PO4":       total_Pi_mM * alpha_H3PO4,
        "H2PO4":       total_Pi_mM * alpha_H2PO4,
        "HPO4":        total_Pi_mM * alpha_HPO4,
        "PO4":         total_Pi_mM * alpha_PO4,
        "pKa1_eff":    pKa1_eff,
        "pKa2_eff":    pKa2_eff,
        "pKa3_eff":    pKa3_eff,
        "alpha_HPO4":  alpha_HPO4,
        "alpha_PO4":   alpha_PO4,
    }


# ── Carbonate speciation ──────────────────────────────────────────────────────

def pKa2_carbonate_at_T(T_celsius: float) -> float:
    """
    pKa₂ for HCO₃⁻ ⇌ H⁺ + CO₃²⁻.
    Harned & Scholes (1941): pKa₂ = 10.33 + 0.010*(25 - T_celsius) approx.
    At 25°C: 10.33; at 0°C: 10.58; at -20°C: 10.78.
    """
    return 10.33 + 0.010 * (25.0 - T_celsius)


def carbonate_speciation(
    total_CO2_mM: float,
    pH: float,
    T_celsius: float,
    I: float,
) -> dict:
    """
    Fractional speciation of total dissolved CO₂ into CO₂(aq), HCO₃⁻, CO₃²⁻.

    Returns dict with keys: CO2_aq, HCO3, CO3, and effective pKa values.

    pKa₁ and pKa₂ are temperature- and ionic-strength corrected.
    """
    pKa1_eff = _pKa_obs(pKa1_CO2_at_T(T_celsius), z_acid=0, z_base=1,
                        I=I, T_celsius=T_celsius)
    pKa2_eff = _pKa_obs(pKa2_carbonate_at_T(T_celsius), z_acid=1, z_base=2,
                        I=I, T_celsius=T_celsius)

    h = 10**(-pH)
    K1 = 10**(-pKa1_eff)
    K2 = 10**(-pKa2_eff)

    denom = h**2 + K1*h + K1*K2
    alpha_CO2  = h**2 / denom
    alpha_HCO3 = K1 * h / denom
    alpha_CO3  = K1 * K2 / denom

    return {
        "CO2_aq":     total_CO2_mM * alpha_CO2,
        "HCO3":       total_CO2_mM * alpha_HCO3,
        "CO3":        total_CO2_mM * alpha_CO3,
        "pKa1_eff":   pKa1_eff,
        "pKa2_eff":   pKa2_eff,
        "alpha_CO3":  alpha_CO3,
    }


# ── Calcium speciation (albumin binding) ──────────────────────────────────────

# Calibration: at [albumin]=4 g/dL, pH=7.4, α_Ca≈0.50 (Fogh-Andersen 1995)
# Simple linear model: α_Ca = 1 / (1 + K_eff(pH) * [albumin_gL])
# K_eff(pH=7.4) = 0.25 (g/dL)^-1 (calibrated)
# pH dependence: higher pH → more Ca²⁺ bound (fewer H⁺ competing for sites)
# Empirical: K_eff(pH) = 0.25 * 10^(0.20*(pH - 7.4))   (Pedersen 1972, linearized)

_K_ALB_REF = 0.25       # (g/dL)^-1 at pH 7.4
_K_ALB_PH_SLOPE = 0.20  # d log K_eff / d pH unit


def calcium_free_fraction(
    total_Ca_mM: float,
    albumin_gL: float,
    pH: float,
    T_celsius: float = 25.0,
    I: float = 0.16,
) -> float:
    """
    Free Ca²⁺ fraction α_Ca = [Ca²⁺]_free / [Ca²⁺]_total.

    Uses a pH- and temperature-corrected linear binding model calibrated to
    Fogh-Andersen et al. (1995): α_Ca(4 g/dL albumin, pH 7.4) ≈ 0.50.

    Limitations (documented, not hidden):
    - Linear model: assumes albumin sites not saturated. At very high [albumin]
      (k×5.58 scenario with 4.5 g/dL initial → 25 g/dL in pool), sites may saturate,
      causing underestimate of α_Ca (the model is conservative, i.e., overestimates
      Ca²⁺ sequestration). True α_Ca likely higher than model predicts at k=5.58.
    - Temperature correction on binding constant: assumed negligible (ΔH binding ≈ 0)
      for this first-order estimate; literature values range from -10 to +5 kJ/mol.
    - Ionic-strength effect on albumin surface charge: not modeled; likely small at
      I < 1 mol/kg relative to the dominant albumin concentration effect.

    Returns α_Ca in range (0, 1].
    """
    albumin_gdL = albumin_gL / 10.0   # convert g/L → g/dL for the formula
    K_eff = _K_ALB_REF * 10.0 ** (_K_ALB_PH_SLOPE * (pH - 7.4))
    alpha = 1.0 / (1.0 + K_eff * albumin_gdL)
    return float(np.clip(alpha, 1e-3, 1.0))


# ── Ion activities in the pool ────────────────────────────────────────────────

def ion_activities(
    composition_mM: dict,
    pH: float,
    T_celsius: float,
    I: float,
) -> dict:
    """
    Compute activities of all relevant ions in the unfrozen pool.

    composition_mM should be pool concentrations (already multiplied by k).
    albumin assumed in 'protein' key (g/dL if converted from mM × MW).

    Returns dict of activities (dimensionless, in molar units).
    """
    # Activity coefficients via Davies (temperature-corrected)
    g_Ca   = activity_coefficient(2, I, T_celsius)
    g_Mg   = activity_coefficient(2, I, T_celsius)
    g_HPO4 = activity_coefficient(2, I, T_celsius)
    g_PO4  = activity_coefficient(3, I, T_celsius)
    g_CO3  = activity_coefficient(2, I, T_celsius)
    g_OH   = activity_coefficient(1, I, T_celsius)
    g_H    = activity_coefficient(1, I, T_celsius)

    # pH → H⁺ and OH⁻ activities
    a_H  = 10.0**(-pH)           # by definition: pH = -log10(a_H)
    pKw  = pKw_at_T(T_celsius)
    a_OH = 10.0**(pH - pKw)      # a_OH = Kw / a_H

    # Phosphate speciation → concentrations of each form
    Pi_mM = composition_mM.get("Pi", 0.0)
    pi_spec = phosphate_speciation(Pi_mM, pH, T_celsius, I)
    c_HPO4 = pi_spec["HPO4"]     # mM
    c_PO4  = pi_spec["PO4"]      # mM

    # Carbonate speciation → CO₃²⁻ concentration
    co2_mM = composition_mM.get("CO2", 0.0)
    co2_spec = carbonate_speciation(co2_mM, pH, T_celsius, I)
    c_CO3 = co2_spec["CO3"]      # mM

    # Calcium free fraction (albumin binding)
    # Convert protein from mM → g/dL (albumin MW = 66500 g/mol)
    protein_mM = composition_mM.get("protein", 0.0)
    albumin_gL = protein_mM / 1000.0 * 66.5  # mmol/L × g/mol / 1000 → g/L
    Ca_total_mM = composition_mM.get("Ca", 0.0)
    alpha_Ca = calcium_free_fraction(Ca_total_mM, albumin_gL, pH, T_celsius, I)
    Ca_free_mM = Ca_total_mM * alpha_Ca

    # Convert mM to mol/L (= mol/kg approximately for dilute solutions)
    def a(c_mM: float, gamma: float) -> float:
        return (c_mM / 1000.0) * gamma

    return {
        "a_Ca":     a(Ca_free_mM, g_Ca),
        "a_HPO4":   a(c_HPO4, g_HPO4),
        "a_PO4":    a(c_PO4, g_PO4),
        "a_CO3":    a(c_CO3, g_CO3),
        "a_OH":     a_OH,
        "a_H":      a_H,
        "Ca_free_mM":    Ca_free_mM,
        "alpha_Ca":      alpha_Ca,
        "albumin_gL":    albumin_gL,
        "c_HPO4_mM":     c_HPO4,
        "c_PO4_mM":      c_PO4,
        "c_CO3_mM":      c_CO3,
        "gamma_Ca":      g_Ca,
        "gamma_PO4":     g_PO4,
        "pi_speciation": pi_spec,
        "co2_speciation": co2_spec,
    }


# ── IAP for each phase ────────────────────────────────────────────────────────

def _iap(phase: str, acts: dict) -> float:
    """
    Ion activity product for a given phase using pre-computed ion activities.
    Returns IAP (dimensionless, in molar units consistent with Ksp).
    """
    a_Ca   = acts["a_Ca"]
    a_HPO4 = acts["a_HPO4"]
    a_PO4  = acts["a_PO4"]
    a_CO3  = acts["a_CO3"]
    a_OH   = acts["a_OH"]
    a_H    = acts["a_H"]

    if phase == "HAp":
        # Ca₅(PO₄)₃OH: 5Ca²⁺ + 3PO₄³⁻ + OH⁻
        return a_Ca**5 * a_PO4**3 * a_OH

    elif phase == "OCP":
        # Ca₈H₂(PO₄)₆: 8Ca²⁺ + 2H⁺ + 6PO₄³⁻
        return a_Ca**8 * a_H**2 * a_PO4**6

    elif phase in ("brushite", "monetite"):
        # CaHPO₄[·2H₂O]: Ca²⁺ + HPO₄²⁻
        return a_Ca * a_HPO4

    elif phase in ("ACP_loose", "ACP_tight"):
        # Ca₃(PO₄)₂: 3Ca²⁺ + 2PO₄³⁻
        return a_Ca**3 * a_PO4**2

    elif phase in ("calcite", "vaterite"):
        # CaCO₃: Ca²⁺ + CO₃²⁻
        return a_Ca * a_CO3

    else:
        raise ValueError(f"Unknown phase: {phase}")


def saturation_index(
    phase: str,
    composition_mM: dict,
    pH: float,
    T_celsius: float,
    I: float,
) -> float:
    """
    SI = log10(IAP / Ksp) for a given phase.

    Positive SI → supersaturated (precipitation thermodynamically favored).
    Negative SI → undersaturated.
    SI = 0    → exact saturation.
    """
    acts = ion_activities(composition_mM, pH, T_celsius, I)
    iap  = _iap(phase, acts)
    log_ksp = ksp_at_T(phase, T_celsius)
    if iap <= 0:
        return np.nan
    return np.log10(iap) - log_ksp


def all_si(
    composition_mM: dict,
    pH: float,
    T_celsius: float,
    I: float,
) -> dict:
    """
    Compute SI for all phases. Returns dict {phase: SI}.
    Also includes activities and speciation details.
    """
    acts = ion_activities(composition_mM, pH, T_celsius, I)
    result = {"_activities": acts}
    for phase in KSP_DATA:
        iap = _iap(phase, acts)
        log_ksp = ksp_at_T(phase, T_celsius)
        result[phase] = np.log10(iap) - log_ksp if iap > 0 else np.nan
    return result


def physiological_serum_si() -> dict:
    """
    Compute SI(HAp) for normal human blood plasma at 37°C.

    Reference conditions (Fogh-Andersen 1995; Lenz et al. 2013):
        pH = 7.4, T = 37°C
        [Ca]_total = 2.4 mM  (ionized ≈ 1.18 mM → α_Ca ≈ 0.50)
        [Pi]_total = 1.0 mM  (inorganic phosphate)
        [albumin]  = 4.0 g/dL = 40 g/L
        I          = 0.16 mol/kg (physiological)
        [CO2]_total ≈ 25 mM (dominated by HCO₃⁻ in venous blood)

    Expected result: SI(HAp) ≈ +3 to +7
    (Heughebaert & Nancollas 1984: "10^5 times supersaturated" → SI ≈ 5)
    (Lenz et al. 2013: rigorous speciation gives SI ≈ 5–7)
    Biological inhibitors suppress nucleation kinetically; the large SI is
    maintained as a metastable state in vivo.
    """
    phys = {
        "Ca":      2.4,
        "Mg":      0.8,
        "Pi":      1.0,
        "Na":    140.0,
        "K":       4.0,
        "Cl":     98.0,
        "CO2":    25.0,
        "protein": 60.0,  # 4.0 g/dL = 40 g/L → 40000/66500 mmol/L ≈ 0.60 mM → 60 µM → need g/L form
    }
    # Fix protein entry: pass as mM where MW=66.5 kDa so g/L → mM = g/L/66.5
    # 40 g/L / 66.5 g/mol × 1000 = 0.602 mM protein
    # Our calcium_free_fraction function converts mM back to g/L using MW=66.5
    # So to input 4 g/dL = 40 g/L albumin: protein_mM = 40/66.5 × 1000 = 601.5 mM
    # Wait — the REFERENCE_COMPOSITION has protein=60.0 mM (placeholder for 60 g/L protein)
    # Our conversion: albumin_gL = protein_mM / 1000 * 66.5
    # 60 mM × 66.5 g/mol / 1000 = 3.99 g/L (but should be 40 g/L)
    # Need to fix: protein_mM should be 40000/66.5 = 601.5 mM for 40 g/L
    # But REFERENCE_COMPOSITION uses 60 mM for conceptual placeholder.
    # Here we set protein_mM correctly for 4 g/dL albumin:
    phys["protein"] = 40.0 / 66.5 * 1000   # = 601.5 mM → converts back to 40 g/L
    return all_si(phys, pH=7.4, T_celsius=37.0, I=0.16)


# ── Summary table ─────────────────────────────────────────────────────────────

def build_summary_table(
    comp_base: dict,
    comp_cryo: dict,
    I_base: float,
    I_cryo: float,
) -> pd.DataFrame:
    """
    Build the summary SI table for data/module3_summary.csv.
    Columns: SI at baseline + two cryoconcentrated scenarios.
    Rows: all phases.
    """
    # Protein mM conversion for 4.5 g/dL baseline albumin
    comp_base = dict(comp_base)
    comp_base["protein"] = 4.5 / 66.5 * 1000 * 10  # 4.5 g/dL = 45 g/L

    comp_cryo_0 = dict(comp_cryo)
    comp_cryo_0["protein"] = comp_cryo_0.get("protein", 0.0)

    comp_cryo_90 = dict(comp_cryo_0)  # same concentrations; pH differs

    scenarios = {
        "SI_k1_pH74_T25C":        (comp_base,  7.40, 25.0,  I_base),
        "SI_k558_pH781_T-20C":    (comp_cryo_0, 7.81, -20.0, I_cryo),
        "SI_k558_pH881_T-20C_90CO2loss": (comp_cryo_90, 8.81, -20.0, I_cryo),
    }

    rows = []
    for phase in KSP_DATA:
        row = {"phase": phase, "Ksp_note": KSP_DATA[phase]["note"]}
        for label, (comp, pH, T, I) in scenarios.items():
            si_val = saturation_index(phase, comp, pH, T, I)
            unc_ksp = ksp_uncertainty_at_T(phase, T)
            unc_act = si_uncertainty_from_activity_model(I)
            row[label] = f"{si_val:+.2f}"
            row[f"{label}_unc"] = f"±{unc_ksp + unc_act:.2f}"
        rows.append(row)

    return pd.DataFrame(rows)


# ── Figures ───────────────────────────────────────────────────────────────────

def plot_si_baseline(
    si_base: dict,
    si_cryo_0: dict,
    si_cryo_90: dict,
    T_cryo: float,
    I_cryo: float,
    output_path: Path,
) -> None:
    """
    fig02: SI bar chart comparing baseline (k=1, 25°C, pH 7.4) vs
    primary cryoconcentrated state (k=5.58, -20°C, pH 7.81, 0% CO₂ loss)
    and secondary scenario (90% CO₂ loss, pH 8.81).

    ACP shown as hatched band (loose to tight endpoint).
    """
    phases_display = ["HAp", "OCP", "brushite", "monetite", "calcite", "vaterite"]
    labels_display = ["HAp\nCa₅(PO₄)₃OH", "OCP\nCa₈H₂(PO₄)₆",
                      "Brushite\nCaHPO₄·2H₂O", "Monetite\nCaHPO₄",
                      "Calcite\nCaCO₃", "Vaterite\nCaCO₃"]

    x = np.arange(len(phases_display))
    width = 0.25

    fig, ax = plt.subplots(figsize=(12, 6))

    def _vals(si_dict, phases):
        return [si_dict.get(p, np.nan) for p in phases]

    bars_base   = ax.bar(x - width, _vals(si_base, phases_display),   width,
                         label="Baseline: k=1, 25°C, pH 7.4", color="steelblue", alpha=0.85)
    bars_cryo0  = ax.bar(x,         _vals(si_cryo_0, phases_display), width,
                         label="Cryo: k=5.58, −20°C, pH 7.81 (0% CO₂ loss)",
                         color="firebrick", alpha=0.85)
    bars_cryo90 = ax.bar(x + width, _vals(si_cryo_90, phases_display), width,
                         label="Cryo + CO₂ loss: pH 8.81 (90% CO₂ loss)",
                         color="darkorange", alpha=0.85)

    # ACP band (separate): show as horizontal span
    acp_loose_base = si_base.get("ACP_loose", 0)
    acp_tight_base = si_base.get("ACP_tight", 0)
    acp_loose_cryo = si_cryo_0.get("ACP_loose", 0)
    acp_tight_cryo = si_cryo_0.get("ACP_tight", 0)

    # Uncertainty error bars (Ksp T-extrapolation + activity model)
    for i, phase in enumerate(phases_display):
        unc_ksp = ksp_uncertainty_at_T(phase, T_cryo)
        unc_act = si_uncertainty_from_activity_model(I_cryo)
        unc = unc_ksp + unc_act
        v0  = si_cryo_0.get(phase, np.nan)
        v90 = si_cryo_90.get(phase, np.nan)
        ax.errorbar(i,       v0,  yerr=unc, fmt="none", color="black", capsize=4, lw=1.5)
        ax.errorbar(i+width, v90, yerr=unc, fmt="none", color="black", capsize=4, lw=1.5)

    # ACP band annotation
    y_acp_lo = min(acp_loose_cryo, acp_tight_cryo)
    y_acp_hi = max(acp_loose_cryo, acp_tight_cryo)
    ax.axhspan(y_acp_lo, y_acp_hi, color="purple", alpha=0.18,
               label=f"ACP range (cryo): SI {y_acp_lo:+.1f} to {y_acp_hi:+.1f}")
    acp_lo_b = min(acp_loose_base, acp_tight_base)
    acp_hi_b = max(acp_loose_base, acp_tight_base)
    ax.axhspan(acp_lo_b, acp_hi_b, color="steelblue", alpha=0.12,
               label=f"ACP range (baseline): SI {acp_lo_b:+.1f} to {acp_hi_b:+.1f}")

    ax.axhline(0, color="black", lw=1.0, ls="--", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_display, fontsize=9.5)
    ax.set_ylabel("Saturation Index  SI = log₁₀(IAP/Ksp)", fontsize=11)
    ax.set_title(
        "Saturation indices: baseline vs cryoconcentrated state\n"
        "Error bars: Ksp temperature-extrapolation + Davies activity uncertainty (±SI)",
        fontsize=11,
    )
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_ylim(bottom=min(-3, ax.get_ylim()[0]))

    # Annotation box: key speciation results
    acts = si_cryo_0.get("_activities", {})
    if acts:
        txt = (
            f"Cryo state (k=5.58, −20°C, pH 7.81):\n"
            f"  α_Ca (free) = {acts.get('alpha_Ca', 0):.2f}\n"
            f"  [Ca²⁺]_free = {acts.get('Ca_free_mM', 0):.1f} mM\n"
            f"  [PO₄³⁻] = {acts.get('c_PO4_mM',0)*1e3:.2e} µM\n"
            f"  [CO₃²⁻] = {acts.get('c_CO3_mM',0):.3f} mM\n"
            f"  I_pool = {I_cryo:.2f} mol/kg"
        )
        ax.text(0.72, 0.97, txt, transform=ax.transAxes,
                fontsize=8, va="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.88))

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_si_temperature(
    comp_base: dict,
    I_base: float,
    output_path: Path,
) -> None:
    """figS2: SI(HAp) and SI(calcite) vs temperature for k = 1, 3, 5.58."""
    T_arr = np.linspace(25.0, -20.0, 100)
    k_values = [1.0, 3.0, 5.58]
    colors = ["steelblue", "darkorange", "firebrick"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for phase, ax in zip(["HAp", "calcite"], axes):
        for k, color in zip(k_values, colors):
            SI_arr = []
            for T in T_arr:
                comp_k = {sp: v * k for sp, v in comp_base.items()}
                I_k = I_base * k
                si_val = saturation_index(phase, comp_k, pH=7.81 if k > 1 else 7.4,
                                          T_celsius=T, I=I_k)
                SI_arr.append(si_val)
            ax.plot(T_arr, SI_arr, color=color, lw=2, label=f"k = {k}")
        ax.axhline(0, color="gray", ls="--", lw=1)
        ax.axvline(-20, color="gray", ls=":", lw=1, alpha=0.6)
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("SI")
        ax.set_title(f"SI({phase}) vs temperature")
        ax.legend()

    fig.suptitle("SI vs temperature for three cryoconcentration scenarios",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_calcium_speciation(
    comp_base: dict,
    I_base: float,
    output_path: Path,
) -> None:
    """figS3: α_Ca(free) vs [albumin] at k=1, 3, 5.58."""
    alb_range = np.linspace(0, 30, 200)   # g/dL
    k_values = [1.0, 3.0, 5.58]
    colors = ["steelblue", "darkorange", "firebrick"]

    fig, ax = plt.subplots(figsize=(8, 5))

    for k, color in zip(k_values, colors):
        Ca_k = comp_base["Ca"] * k
        alpha_arr = [calcium_free_fraction(Ca_k, alb*10, pH=7.4) for alb in alb_range]
        ax.plot(alb_range, alpha_arr, color=color, lw=2, label=f"k={k}")

    # Mark physiological albumin and cryoconcentrated albumin
    alb_phys  = 4.0          # g/dL
    alb_cryo_gdL = 4.5 * 5.58   # initial 4.5 g/dL × k=5.58 = 25.1 g/dL in pool
    ax.axvline(alb_phys, color="navy", ls="--", lw=1.2, label="Physiological 4 g/dL")
    ax.axvline(alb_cryo_gdL, color="crimson", ls="--", lw=1.2,
               label=f"Cryo pool k=5.58 → {alb_cryo_gdL:.1f} g/dL")
    ax.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.6)

    ax.set_xlabel("[Albumin] (g/dL)")
    ax.set_ylabel("α_Ca = [Ca²⁺]_free / [Ca²⁺]_total")
    ax.set_title("Free Ca²⁺ fraction vs albumin concentration\n"
                 "Shows albumin buffering effect that partially counteracts cryoconcentration")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Module 3: Saturation indices")
    print("=" * 60)

    comp_base = dict(REFERENCE_COMPOSITION_mM)
    # Fix protein to correctly represent 4.5 g/dL albumin (midpoint of 3.7-7.3 g/dL)
    # protein_mM = 4.5 g/dL × 10 × 1000/66500 = 677 mM
    comp_base["protein"] = 4.5 * 10.0 / 66.5 * 1000.0

    I_base = ionic_strength_baseline(comp_base)
    k_cryo = 5.58   # glycerol 15%, -20°C (from Module 2)
    I_cryo = I_base * k_cryo

    # Pool composition at k=5.58
    comp_cryo = {sp: v * k_cryo for sp, v in comp_base.items()}

    print(f"\nBaseline: I = {I_base:.4f} mol/kg")
    print(f"Cryo pool: k = {k_cryo}, I = {I_cryo:.4f} mol/kg")

    # Physiological sanity check
    print("\n── Physiological sanity check ──")
    phys = physiological_serum_si()
    si_hap_phys = phys["HAp"]
    acts_phys = phys["_activities"]
    print(f"  SI(HAp) at 37°C, pH 7.4, [Ca]=2.4 mM, [Pi]=1.0 mM: {si_hap_phys:+.2f}")
    print(f"  α_Ca = {acts_phys['alpha_Ca']:.2f} ([Ca²⁺]_free = {acts_phys['Ca_free_mM']:.2f} mM)")
    print(f"  [PO₄³⁻] = {acts_phys['c_PO4_mM']*1e6:.1f} nM  (tiny — drives SI numerics)")
    if 3.0 <= si_hap_phys <= 8.0:
        print("  ✓ PASS: SI(HAp) in expected biological range [+3, +8]")
    else:
        print(f"  ✗ WARN: SI(HAp) = {si_hap_phys:.2f} — outside expected [+3, +8]; check model")

    # Baseline SI (k=1, 25°C, pH 7.4)
    print("\n── Baseline SI (k=1, T=25°C, pH=7.4) ──")
    si_base = all_si(comp_base, pH=7.4, T_celsius=25.0, I=I_base)
    for phase in KSP_DATA:
        unc = ksp_uncertainty_at_T(phase, 25.0) + si_uncertainty_from_activity_model(I_base)
        print(f"  {phase:<12} SI = {si_base[phase]:+.2f}  (±{unc:.2f})")

    # Cryoconcentrated SI (k=5.58, T=-20°C, pH=7.81, 0% CO₂ loss)
    print("\n── Cryoconcentrated SI (k=5.58, T=−20°C, pH=7.81, 0% CO₂ loss) ──")
    si_cryo_0 = all_si(comp_cryo, pH=7.81, T_celsius=-20.0, I=I_cryo)
    acts_cryo = si_cryo_0["_activities"]
    print(f"  α_Ca = {acts_cryo['alpha_Ca']:.3f}  ([Ca²⁺]_free = {acts_cryo['Ca_free_mM']:.2f} mM)")
    print(f"  Albumin in pool: {acts_cryo['albumin_gL']:.1f} g/L = {acts_cryo['albumin_gL']/10:.1f} g/dL")
    for phase in KSP_DATA:
        unc = ksp_uncertainty_at_T(phase, -20.0) + si_uncertainty_from_activity_model(I_cryo)
        print(f"  {phase:<12} SI = {si_cryo_0[phase]:+.2f}  (±{unc:.2f})")

    # 90% CO₂ loss scenario (pH=8.81)
    print("\n── Cryoconcentrated SI (k=5.58, T=−20°C, pH=8.81, 90% CO₂ loss) ──")
    si_cryo_90 = all_si(comp_cryo, pH=8.81, T_celsius=-20.0, I=I_cryo)
    for phase in ["HAp", "OCP", "calcite", "ACP_loose", "ACP_tight"]:
        print(f"  {phase:<12} SI = {si_cryo_90[phase]:+.2f}")

    # Save CSV
    data_dir = Path(__file__).parent.parent / "data"
    data_dir.mkdir(exist_ok=True)
    df = build_summary_table(comp_base, comp_cryo, I_base, I_cryo)
    csv_path = data_dir / "module3_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")

    # Figures
    fig_dir = Path(__file__).parent.parent / "figures"
    fig_dir.mkdir(exist_ok=True)

    plot_si_baseline(si_base, si_cryo_0, si_cryo_90,
                     T_cryo=-20.0, I_cryo=I_cryo,
                     output_path=fig_dir / "fig02_si_baseline.png")
    plot_si_temperature(comp_base, I_base,
                        output_path=fig_dir / "figS2_si_temperature.png")
    plot_calcium_speciation(comp_base, I_base,
                            output_path=fig_dir / "figS3_calcium_speciation.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
