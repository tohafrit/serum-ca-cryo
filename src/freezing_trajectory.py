"""
Module 2 — Cryoconcentration trajectory and pH evolution in frozen serum QC standard.

Physical story
--------------
When serum freezes, pure ice excludes solutes. All ions concentrate in a shrinking
unfrozen interstitial pool. The concentration multiplier k = 1/f, where f is the
unfrozen water fraction. k is bounded by the initial freezing-point depression:

    f(T) = ΔTf_initial / |T|   (ideal dilute colligative model)

Even at the conservative k ≈ 5 (glycerol-protected matrix), the Ca·Pi ion product
increases ~23-fold — sufficient to cross supersaturation thresholds for all calcium
phosphate phases (verified in Module 3).

pH evolution: pure cryoconcentration does NOT change pH (both CO₂ and HCO₃⁻ scale
equally). pH rises come from (a) CO₂ outgassing to headspace/ice, (b) increasing
pKa₁(CO₂) at low temperature, (c) ionic strength effects on effective pKa.

Validity of ideal van't Hoff
-----------------------------
The ideal dilute approximation breaks down when ionic strength I > ~0.5 mol/kg
(not at total osmolality, because glycerol is a non-electrolyte and does not
disturb ionic activity coefficients nearly as strongly as salts at the same
osmolality). Threshold:

    I < 0.5 mol/kg  → Davies equation valid (OK)
    0.5–2.0         → Extended D-H / Davies approximate (±20–30%)
    > 2.0           → Pitzer model required (pyEQL)

References
----------
- Kf for water: Atkins & de Paula, Physical Chemistry 10th ed., p. 151
- CO₂ pKa vs temperature: Harned & Davis (1943) J Am Chem Soc 65:2030–2037
- Phosphate pKa vs temperature: Goldberg et al. (2002) J Phys Chem Ref Data 31:231–370
- Davies activity correction: Davies CW (1962) Ion Association, Butterworths
- Glycerol osmolality and cryoprotection: Pegg DE (2002) Cryobiology 44:58–69
- Glycerol viscosity at low T: Segur & Oberstar (1951) Ind Eng Chem 43:2117–2120
- DMSO-water thermodynamics: Cowie & Toporowski (1961) Can J Chem 39:2240–2243
"""

from __future__ import annotations

from typing import Literal
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Physical constants ────────────────────────────────────────────────────────

Kf_WATER = 1.853  # cryoscopic constant K·kg/mol (Atkins 10th ed.)
A_DAVIES_25C = 0.509  # Davies A coefficient at 25°C (dimensionless)

# van't Hoff factors for osmolality (all species treated as fully dissociated)
# protein VHF=0 because albumin (~0.9 mM) adds <1 mOsm/kg — negligible
VHF = {
    "Na": 1, "K": 1, "Cl": 1, "Ca": 1, "Mg": 1,
    "Pi": 1, "CO2": 2, "protein": 0,
}

# Effective charge² for ionic strength: Pi and CO2 computed with pH-dependent
# speciation; others use integer charges
ION_CHARGE_SQ = {
    "Na": 1, "K": 1, "Cl": 1,
    "Ca": 4,   # z=2 → z²=4
    "Mg": 4,
    # Pi and CO2: pH-dependent, handled in ionic_strength_baseline()
    "protein": 0,
}

# Cryoprotectant osmolality contributions (mmol/kg = mOsm/kg, non-electrolyte i=1)
# Glycerol 15% w/w: 150 g/kg / 92.09 g/mol × 1000 = 1629 mmol/kg; rounded to 1925
# to partially account for non-ideal solute behavior at this concentration
# DMSO 10% w/w: 100 g/kg / 78.13 g/mol × 1000 = 1280 mmol/kg
CRYO_OSMOLALITY_mOsm_PER_KG = {
    "glycerol_15pct": 1629.0,
    "dmso_10pct":     1280.0,
    "none":              0.0,
}

CRYO_LABEL = {
    "glycerol_15pct": "Glycerol 15% w/w",
    "dmso_10pct":     "DMSO 10% w/w",
    "none":           "No cryoprotectant",
}

# Viscosity of unfrozen phase at −20°C relative to water at 20°C (order of magnitude)
# Used only for kinetic arguments in Module 5; noted here for cross-module consistency
CRYO_VISCOSITY_FACTOR_MINUS20 = {
    "glycerol_15pct": 500.0,  # Segur & Oberstar (1951)
    "dmso_10pct":      50.0,  # estimate from DMSO-water data
    "none":            20.0,  # supercooled pure water
}

# Reference serum formulation — mid-range of challenge brief values
REFERENCE_COMPOSITION_mM = {
    "Na":      140.0,   # 118–160 mM
    "K":         5.0,   # 3.5–7.2 mM
    "Cl":       93.0,   # 82–104 mM
    "Ca":        2.75,  # 6.6–15.6 mg/dL ÷ 40.08 g/mol
    "Mg":        1.25,  # 1.7–4.6 mg/dL ÷ 24.31 g/mol
    "Pi":        1.80,  # 2.3–8.8 mg/dL as P ÷ 30.97
    "CO2":      30.0,   # 16–45 mmol/L
    "protein":  60.0,   # ~60 g/L; included for completeness, VHF=0
}


# ── Temperature-dependent pKa values ─────────────────────────────────────────

def pKa1_CO2_at_T(T_celsius: float) -> float:
    """
    Apparent pKa₁ for CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻ at temperature T.

    Linear fit to Harned & Davis (1943) data (valid 0–40°C; extrapolated to −20°C):
        pKa₁ = 6.352 + 0.009 × (25 − T°C)
    At 25°C: 6.352; at 0°C: 6.577 (lit. 6.581); at 37°C: 6.244 (lit. 6.284).
    Extrapolation to −20°C: 6.757 (order-of-magnitude confidence).
    """
    return 6.352 + 0.009 * (25.0 - T_celsius)


def pKa2_phosphate_at_T(T_celsius: float) -> float:
    """
    Apparent pKa₂ for H₂PO₄⁻ ⇌ H⁺ + HPO₄²⁻ at temperature T.

    Linear fit to Goldberg et al. (2002) data:
        pKa₂ = 7.21 + 0.0115 × (25 − T°C)
    At 25°C: 7.21; at 0°C: 7.50 (lit. 7.54); at 37°C: 7.07 (lit. 7.07).
    Extrapolation to −20°C: 7.73 (order-of-magnitude confidence).
    """
    return 7.21 + 0.0115 * (25.0 - T_celsius)


# ── Davies activity correction ────────────────────────────────────────────────

def davies_log_gamma(z: float, I: float, A: float = A_DAVIES_25C) -> float:
    """
    log₁₀(γ) for an ion of charge z at ionic strength I (mol/kg) via Davies equation.

    Davies (1962): log₁₀(γ) = −A·z²·(√I/(1+√I) − 0.3·I)
    Valid for I ≲ 0.5 mol/kg; increasingly approximate above that.
    A ≈ 0.509 at 25°C; temperature dependence small over 0–25°C range.

    Returns a negative number (γ < 1 for all ions at I > 0).
    """
    sqrtI = np.sqrt(max(I, 1e-12))
    return -A * z**2 * (sqrtI / (1.0 + sqrtI) - 0.3 * I)


def pKa_effective(pKa_thermo: float, z_product: float, I: float) -> float:
    """
    Effective (observed) pKa corrected for ionic strength via Davies equation.

    For reaction HA → H⁺ + A^(z_product):
        pKa_eff = pKa_thermo − log₁₀(γ_H+ · γ_A / γ_HA)
        γ_HA = 1 if HA is neutral (CO₂, H₂O)
        For HA charged: adjust accordingly.

    Parameters
    ----------
    pKa_thermo : float  — thermodynamic pKa (dilute limit)
    z_product  : float  — charge of the anion product A^z
    I          : float  — ionic strength (mol/kg)
    """
    log_g_H = davies_log_gamma(1, I)
    log_g_A = davies_log_gamma(abs(z_product), I)
    # For CO₂/H₃PO₄ as HA (neutral or singly charged), γ_HA correction:
    # Neutral HA: log_g_HA = 0 → correction = log_g_H + log_g_A
    # Singly negative HA (H₂PO₄⁻ → HPO₄²⁻): log_g_HA = log_g_z1
    if abs(z_product) == 1:
        # HA is neutral (CO₂ case), product is monovalent
        log_g_HA = 0.0
    else:
        # HA is monovalent (H₂PO₄⁻ case), product is divalent
        log_g_HA = davies_log_gamma(1, I)
    return pKa_thermo - (log_g_H + log_g_A - log_g_HA)


# ── Ionic strength ────────────────────────────────────────────────────────────

def ionic_strength_baseline(
    composition_mM: dict,
    pH: float = 7.4,
) -> float:
    """
    Baseline ionic strength (mol/kg) at k=1 using I = 0.5·Σ(cᵢ·zᵢ²).

    Pi and CO₂ contribute pH-dependent fractions of their species.
    Uses pKa values at 25°C for baseline.

    Result ~0.14–0.16 mol/kg for typical serum (consistent with known ~0.15 mol/L).
    """
    comp = {k: v / 1000.0 for k, v in composition_mM.items()}  # mmol/L → mol/L ≈ mol/kg

    I = 0.0
    for ion, z2 in ION_CHARGE_SQ.items():
        if ion in comp:
            I += 0.5 * comp[ion] * z2

    # Pi: mix of H₂PO₄⁻ (z=1) and HPO₄²⁻ (z=2) at given pH
    if "Pi" in comp:
        pKa2 = 7.21  # at 25°C
        f_HPO4 = 10**(pH - pKa2) / (1 + 10**(pH - pKa2))
        f_H2PO4 = 1 - f_HPO4
        z2_Pi = f_H2PO4 * 1 + f_HPO4 * 4  # mean z² weighted by fraction
        I += 0.5 * comp["Pi"] * z2_Pi

    # CO₂: mostly HCO₃⁻ (z=1) at pH 7.4; CO₃²⁻ negligible (pKa₂=10.33)
    if "CO2" in comp:
        pKa1 = 6.352  # at 25°C
        f_HCO3 = 10**(pH - pKa1) / (1 + 10**(pH - pKa1))
        I += 0.5 * comp["CO2"] * f_HCO3 * 1  # z²=1 for HCO₃⁻

    return I


# ── Validity checker ──────────────────────────────────────────────────────────

def assumption_validity_check(k: float, ionic_str_baseline_mol_per_kg: float) -> dict:
    """
    Rate validity of ideal van't Hoff based on ionic strength of the pool.

    Ionic strength (not total osmolality) governs activity coefficient deviations
    for electrolytes. Glycerol is a non-electrolyte: its high osmolality does NOT
    make ionic thermodynamics non-ideal (though it slightly reduces ε of water,
    affecting the Davies A coefficient by <10% at 15% glycerol).

    Zones (empirical, Pitzer 1973 + Davies 1962):
        I < 0.5 mol/kg  → Davies equation valid (<5% error in γ±)    OK
        0.5–2.0         → Extended D-H / Davies approximate           approximate
        > 2.0           → Pitzer model needed for accuracy            unreliable
    """
    I_pool = ionic_str_baseline_mol_per_kg * k
    if I_pool < 0.5:
        level, color = "OK", "green"
        msg = f"I={I_pool:.2f} mol/kg — Davies valid (<5% error in γ±)."
    elif I_pool < 2.0:
        level, color = "approximate", "orange"
        msg = f"I={I_pool:.2f} mol/kg — approximate (±20–30%); Pitzer recommended."
    else:
        level, color = "unreliable", "red"
        msg = f"I={I_pool:.2f} mol/kg — strongly non-ideal; use Pitzer (pyEQL)."
    return {"level": level, "color": color, "I_pool": I_pool, "message": msg}


# ── Core cryoconcentration functions ─────────────────────────────────────────

def freezing_point_depression(
    composition_mM: dict,
    cryoprotectant: Literal["glycerol_15pct", "dmso_10pct", "none"] = "glycerol_15pct",
) -> dict:
    """
    Freezing-point depression ΔTf at k=1 (initial formulation).

    ΔTf = Kf · b_eff   where b_eff = total osmolality (mol/kg).
    Glycerol/DMSO are non-electrolytes (i=1); they contribute osmolality
    but not ionic strength.

    Returns
    -------
    dict: delta_T, osm_electrolyte, osm_cryo, osm_total (all Osm/kg or °C)
    """
    osm_elec = sum(
        composition_mM.get(k, 0.0) * VHF.get(k, 1) / 1000.0
        for k in VHF
    )
    osm_cryo = CRYO_OSMOLALITY_mOsm_PER_KG[cryoprotectant] / 1000.0
    osm_total = osm_elec + osm_cryo
    return {
        "delta_T": Kf_WATER * osm_total,
        "osm_electrolyte": osm_elec,
        "osm_cryo": osm_cryo,
        "osm_total": osm_total,
        "cryoprotectant": cryoprotectant,
    }


def unfrozen_fraction(
    T_celsius: float | np.ndarray,
    fpd_result: dict,
) -> np.ndarray:
    """
    Fraction of initial water remaining liquid at temperature T ≤ 0°C.

    Derivation: at equilibrium, ΔTf(pool) = ΔTf_initial / f → f = ΔTf_initial / |T|.
    Clipped to [0, 1]. At T ≥ 0: f = 1 (all liquid).
    """
    T = np.atleast_1d(np.asarray(T_celsius, dtype=float))
    delta_T0 = fpd_result["delta_T"]
    T_safe = np.where(T == 0.0, -1e-30, T)
    f = np.where(
        T >= 0,
        1.0,
        np.clip(delta_T0 / np.abs(T_safe), 0.0, 1.0),
    )
    return f


def cryoconcentration_trajectory(
    T_array: np.ndarray,
    composition_mM: dict,
    cryoprotectant: Literal["glycerol_15pct", "dmso_10pct", "none"] = "glycerol_15pct",
) -> pd.DataFrame:
    """
    Cryoconcentration trajectory from 0°C to min(T_array).

    Returns a DataFrame with T, f, k, ionic strength, validity level,
    and concentrations of all species in the unfrozen pool.

    Notes
    -----
    k is capped at 100 to avoid Inf near eutectic.
    Ionic strength at each k = I_baseline × k (linear; valid since all species
    concentrate equally in the ideal model).
    """
    fpd = freezing_point_depression(composition_mM, cryoprotectant)
    I_base = ionic_strength_baseline(composition_mM)
    T_arr = np.asarray(T_array, dtype=float)
    f_arr = unfrozen_fraction(T_arr, fpd)
    k_arr = np.clip(1.0 / np.where(f_arr > 0, f_arr, 1e-9), 1.0, 100.0)

    rows = []
    for T, f, k in zip(T_arr, f_arr, k_arr):
        v = assumption_validity_check(k, I_base)
        row = {
            "T_celsius": T,
            "unfrozen_fraction": f,
            "k": k,
            "ionic_strength_pool": I_base * k,
            "validity_level": v["level"],
        }
        for species, c0 in composition_mM.items():
            row[f"{species}_mM"] = c0 * k
        rows.append(row)

    return pd.DataFrame(rows)


# ── pH trajectory ─────────────────────────────────────────────────────────────

def pH_trajectory(
    k_array: np.ndarray,
    composition_mM: dict,
    co2_loss_fractions: tuple = (0.0, 0.5, 0.9),
    T_celsius: float = -20.0,
    pH_initial: float = 7.4,
) -> pd.DataFrame:
    """
    Predicted pH in unfrozen pool vs cryoconcentration factor k,
    for different CO₂ loss fractions during freezing.

    Physical basis
    --------------
    Pure cryoconcentration does NOT change pH: both CO₂ and HCO₃⁻ scale by k,
    the ratio is unchanged, and H-H gives the same pH. pH shifts come from:

    1. Temperature: pKa₁(CO₂) increases at lower T (+0.41 units at −20°C vs 25°C)
       → at same [HCO₃⁻]/[CO₂] ratio, pH is higher.
    2. CO₂ outgassing: CO₂ escapes to headspace/ice more than HCO₃⁻ →
       ratio [HCO₃⁻]/[CO₂] rises → pH rises. Each halving of CO₂ raises pH by +0.30.
    3. Davies ionic strength correction to pKa₁ (secondary, ~+0.25 units).

    H-H equation: pH = pKa₁_eff(T, I) + log₁₀([HCO₃⁻]_pool / [CO₂]_pool)

    At initial conditions (k=1, T=25°C), pH_initial defines:
        log₁₀([HCO₃⁻]₀/[CO₂]₀) = pH_initial − pKa₁_eff(25°C, I_baseline)

    At state (k, T, f_loss):
        [CO₂]_pool = k·[CO₂]₀·(1−f_loss)
        [HCO₃⁻]_pool = k·[HCO₃⁻]₀
        pH = pKa₁_eff(T, k·I_base) + log₁₀_ratio_initial − log₁₀(1−f_loss)

    Note: k cancels in the ratio (but enters through ionic strength correction).

    Parameters
    ----------
    k_array          : concentration factors to evaluate
    composition_mM   : initial composition; needs 'CO2' key
    co2_loss_fractions : fractions of CO₂(aq) lost to headspace/ice (0=none, 1=all)
    T_celsius        : temperature of unfrozen pool (default −20°C)
    pH_initial       : pH of formulation at 25°C, k=1

    Returns
    -------
    DataFrame with columns: k, and one pH column per co2_loss_fraction.
    Also includes pKa1_eff and ionic_strength_pool columns.
    """
    I_base = ionic_strength_baseline(composition_mM, pH=pH_initial)
    k_arr = np.asarray(k_array, dtype=float)

    rows = []
    for k in k_arr:
        I_pool = I_base * k
        pKa1_T = pKa1_CO2_at_T(T_celsius)
        pKa1_eff = pKa_effective(pKa1_T, z_product=1, I=I_pool)

        # Baseline ratio defined by initial pH at 25°C with ionic strength correction
        pKa1_eff_25 = pKa_effective(pKa1_CO2_at_T(25.0), z_product=1, I=I_base)
        log_ratio_initial = pH_initial - pKa1_eff_25

        row = {"k": k, "pKa1_eff": pKa1_eff, "ionic_strength_pool": I_pool}
        for f_loss in co2_loss_fractions:
            if f_loss >= 1.0:
                pH = np.nan  # all CO₂ gone → pH undefined by H-H
            else:
                log_ratio = log_ratio_initial - np.log10(1.0 - f_loss)
                pH = pKa1_eff + log_ratio
            col = f"pH_co2loss_{int(f_loss*100):02d}pct"
            row[col] = pH
        rows.append(row)

    return pd.DataFrame(rows)


# ── Figure: fig01 multi-scenario cryoconcentration ────────────────────────────

def plot_cryoconcentration(
    scenarios: dict,
    primary_key: str,
    I_base: float,
    output_path: Path,
) -> None:
    """
    Dual-axis plot of unfrozen fraction vs temperature for multiple cryoprotectant
    scenarios. Background shading based on ionic strength validity of primary scenario.

    scenarios: dict of {cryo_key: (df, fpd)}
    """
    COLORS = {
        "glycerol_15pct": ("steelblue",  "firebrick"),
        "dmso_10pct":     ("darkorange", "saddlebrown"),
        "none":           ("purple",     "darkmagenta"),
    }
    STYLES = {
        "glycerol_15pct": "-",
        "dmso_10pct":     "--",
        "none":           ":",
    }

    fig, ax1 = plt.subplots(figsize=(10, 5.8))
    ax2 = ax1.twinx()

    # Shading based on primary scenario ionic strength
    primary_df, primary_fpd = scenarios[primary_key]
    T_arr = primary_df["T_celsius"].values

    # Find validity transition temperatures (ionic strength based)
    v_levels = primary_df["validity_level"].values
    t_approx = t_unrel = None
    for i, v in enumerate(v_levels):
        if v == "approximate" and t_approx is None:
            t_approx = T_arr[i]
        if v == "unreliable" and t_unrel is None:
            t_unrel = T_arr[i]

    t_start = T_arr[0]
    t_end = T_arr[-1]
    ax1.axvspan(t_start, t_approx or t_end, color="#c8e6c9", alpha=0.35, zorder=0)
    if t_approx:
        ax1.axvspan(t_approx, t_unrel or t_end, color="#fff9c4", alpha=0.5, zorder=0)
    if t_unrel:
        ax1.axvspan(t_unrel, t_end, color="#ffcdd2", alpha=0.5, zorder=0)

    # Plot each scenario
    for cryo_key, (df, fpd) in scenarios.items():
        T = df["T_celsius"].values
        f = df["unfrozen_fraction"].values
        k = df["k"].values
        c_f, c_k = COLORS.get(cryo_key, ("gray", "gray"))
        ls = STYLES.get(cryo_key, "-")
        label = CRYO_LABEL.get(cryo_key, cryo_key)
        k_final = k[-1]

        ax1.plot(T, f, color=c_f, lw=2.2, ls=ls, label=f"$f$ — {label}")
        ax2.plot(T, k, color=c_k, lw=2.2, ls=ls, label=f"$k$ — {label}")

        # Annotate k at −20°C
        ax2.annotate(
            f"$k$={k_final:.1f}",
            xy=(T[-1], k_final),
            xytext=(T[-1] - 2.5, k_final + 0.3),
            color=c_k, fontsize=9, ha="right",
            arrowprops=dict(arrowstyle="->", color=c_k, lw=0.8),
        )

    # Axes labels
    ax1.set_xlabel("Temperature (°C)", fontsize=12)
    ax1.set_ylabel("Unfrozen water fraction $f$", color="steelblue", fontsize=11)
    ax1.tick_params(axis="y", labelcolor="steelblue")
    ax1.set_xlim(t_start, t_end)
    ax1.set_ylim(-0.05, 1.15)
    ax2.set_ylabel("Concentration multiplier $k = 1/f$", color="firebrick", fontsize=11)
    ax2.tick_params(axis="y", labelcolor="firebrick")

    # Validity legend patches
    ok_patch    = mpatches.Patch(color="#c8e6c9", alpha=0.7, label="OK (Davies valid, I<0.5)")
    apx_patch   = mpatches.Patch(color="#fff9c4", alpha=0.7, label="Approximate (0.5<I<2.0)")
    bad_patch   = mpatches.Patch(color="#ffcdd2", alpha=0.7, label="Unreliable (I>2.0, Pitzer)")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(
        h1 + h2 + [ok_patch, apx_patch, bad_patch],
        l1 + l2 + ["OK (Davies)", "Approximate", "Unreliable (Pitzer)"],
        loc="upper right", fontsize=8.5, ncol=2,
    )

    # Info box: Ca·Pi product at k_final for primary scenario
    row_final = scenarios[primary_key][0].iloc[-1]
    k_prim = row_final["k"]
    ca0, pi0 = REFERENCE_COMPOSITION_mM["Ca"], REFERENCE_COMPOSITION_mM["Pi"]
    iap_ratio = (ca0 * k_prim) * (pi0 * k_prim) / (ca0 * pi0)
    textstr = (
        f"Primary scenario (glycerol 15%):\n"
        f"  ΔTf = {scenarios[primary_key][1]['delta_T']:.2f} °C\n"
        f"  k at −20°C = {k_prim:.2f}×\n"
        f"  Ca·Pi product × {iap_ratio:.0f} (vs k=1)\n"
        f"  I_pool = {row_final['ionic_strength_pool']:.2f} mol/kg"
    )
    ax1.text(0.02, 0.42, textstr, transform=ax1.transAxes,
             fontsize=8.2, va="top",
             bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    ax1.set_title(
        "Cryoconcentration trajectory — serum QC standard\n"
        "Three cryoprotectant scenarios; shading = ionic strength validity of ideal van't Hoff",
        fontsize=11,
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ── Figure: figS1 pH trajectory ───────────────────────────────────────────────

def plot_pH_trajectory(
    df_pH: pd.DataFrame,
    k_markers: dict,
    output_path: Path,
) -> None:
    """
    figS1: pH in unfrozen pool vs k for three CO₂ loss scenarios.
    k_markers: dict of {label: k_value} to mark cryoprotectant scenarios.
    """
    CO2_LOSS_COLS = {
        "pH_co2loss_00pct": ("0% CO₂ loss (sealed vial)", "steelblue",  "-"),
        "pH_co2loss_50pct": ("50% CO₂ loss",              "darkorange", "--"),
        "pH_co2loss_90pct": ("90% CO₂ loss",              "crimson",    ":"),
    }

    fig, ax = plt.subplots(figsize=(9, 5.5))

    for col, (label, color, ls) in CO2_LOSS_COLS.items():
        if col in df_pH.columns:
            ax.plot(df_pH["k"], df_pH[col], color=color, lw=2.5, ls=ls, label=label)

    # Mark cryoprotectant k values
    for label, k_val in k_markers.items():
        ax.axvline(x=k_val, color="gray", ls="-.", lw=1.2, alpha=0.7)
        ax.text(k_val + 0.4, 7.45, label, fontsize=8.5, color="gray",
                rotation=90, va="bottom")

    # Shade CaP supersaturation regions
    ax.axhspan(7.4, 8.0,  color="#fff9c4", alpha=0.5, label="Moderate CaP supersaturation")
    ax.axhspan(8.0, 10.0, color="#ffcdd2", alpha=0.4, label="High CaP supersaturation")
    ax.axhline(y=7.4, color="navy", ls="--", lw=1.2, alpha=0.7, label="Initial pH 7.4")

    ax.set_xlabel("Cryoconcentration factor $k$", fontsize=12)
    ax.set_ylabel("pH in unfrozen pool", fontsize=12)
    ax.set_xlim(1, df_pH["k"].max())
    ax.set_ylim(7.3, 9.2)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title(
        "pH trajectory in unfrozen pool during freezing\n"
        f"T = −20°C; accounting for temperature correction of pKa₁(CO₂) and Davies ionic strength",
        fontsize=11,
    )

    # Annotation: dominant effects
    ax.text(0.62, 0.18,
            "pH rise sources:\n"
            "①  pKa₁(CO₂): +0.41 at −20°C vs 25°C\n"
            "②  CO₂ loss: +0.30 per halving of [CO₂]\n"
            "③  Davies I correction: +0.25 (secondary)",
            transform=ax.transAxes, fontsize=8.5,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Module 2: Cryoconcentration trajectory (multi-scenario + pH)")
    print("=" * 60)

    comp = REFERENCE_COMPOSITION_mM
    I_base = ionic_strength_baseline(comp)
    print(f"\nBaseline ionic strength: {I_base:.4f} mol/kg")

    # ── Scenario trajectories ─────────────────────────────────────────────────
    T_arr = np.linspace(-0.01, -20.0, 500)
    cryo_keys = ["glycerol_15pct", "dmso_10pct", "none"]
    scenarios = {}

    print(f"\n{'Scenario':<22} {'ΔTf (°C)':<12} {'k at -20°C':<14} {'I_pool (mol/kg)'}")
    print("-" * 65)
    for cryo in cryo_keys:
        fpd = freezing_point_depression(comp, cryo)
        df = cryoconcentration_trajectory(T_arr, comp, cryo)
        scenarios[cryo] = (df, fpd)
        row_20 = df.iloc[-1]
        print(f"  {CRYO_LABEL[cryo]:<20} {fpd['delta_T']:<12.3f} "
              f"{row_20['k']:<14.2f} {row_20['ionic_strength_pool']:.3f}")

    # Validity check for primary scenario at k_final
    primary_df, primary_fpd = scenarios["glycerol_15pct"]
    k_glycerol = primary_df.iloc[-1]["k"]
    v = assumption_validity_check(k_glycerol, I_base)
    print(f"\nValidity at k={k_glycerol:.1f} (glycerol 15%): {v['level'].upper()}")
    print(f"  → {v['message']}")

    # Ca·Pi product scaling
    ca0, pi0 = comp["Ca"], comp["Pi"]
    iap_ratio = (k_glycerol**2)
    print(f"\nCa·Pi ion product at k={k_glycerol:.1f}: ×{iap_ratio:.1f} vs baseline")
    print(f"  Ca in pool: {ca0*k_glycerol:.1f} mM (from {ca0} mM)")
    print(f"  Pi in pool: {pi0*k_glycerol:.1f} mM (from {pi0} mM)")

    # ── pH trajectory ─────────────────────────────────────────────────────────
    k_for_pH = np.linspace(1.0, 40.0, 300)
    df_pH = pH_trajectory(k_for_pH, comp, co2_loss_fractions=(0.0, 0.5, 0.9),
                          T_celsius=-20.0)
    k_dmso = scenarios["dmso_10pct"][0].iloc[-1]["k"]
    k_none = scenarios["none"][0].iloc[-1]["k"]

    print(f"\npH at k={k_glycerol:.1f} (glycerol 15%), T=-20°C:")
    row_pH = df_pH.iloc[(df_pH["k"] - k_glycerol).abs().argmin()]
    print(f"  0%  CO₂ loss: {row_pH['pH_co2loss_00pct']:.2f}")
    print(f"  50% CO₂ loss: {row_pH['pH_co2loss_50pct']:.2f}")
    print(f"  90% CO₂ loss: {row_pH['pH_co2loss_90pct']:.2f}")

    # ── Generate figures ──────────────────────────────────────────────────────
    fig_dir = Path(__file__).parent.parent / "figures"
    fig_dir.mkdir(exist_ok=True)

    plot_cryoconcentration(
        scenarios, "glycerol_15pct", I_base,
        fig_dir / "fig01_cryoconcentration.png",
    )
    plot_pH_trajectory(
        df_pH,
        k_markers={
            f"Glycerol 15%\nk={k_glycerol:.1f}": k_glycerol,
            f"DMSO 10%\nk={k_dmso:.1f}": k_dmso,
            f"No cryo\nk={k_none:.1f}": k_none,
        },
        output_path=fig_dir / "figS1_ph_trajectory.png",
    )
    print("\nDone.")


if __name__ == "__main__":
    main()
