"""
Module 7: Intervention modeling.

Three formulation-neutral interventions evaluated against a realistic
baseline (pH=8.0 at −20°C, representing 30% CO₂ outgassing typical in
practice).  Module 6 used the optimistic pH=7.81 (0% loss); Module 7
establishes the true baseline and quantifies each intervention's effect.

Interventions:
  A — Pre-freeze vacuum degassing   (prevents CO₂-driven pH rise)
  B — Controlled-rate freezing      (reduces k heterogeneity + delays nucleation)
  C — Modified thaw protocol        (vortex reduces Noyes-Whitney boundary layer)
  Combined — A + B + C together

Compliance:  all interventions are formulation-neutral.
  A: removes dissolved gas only; chemical composition unchanged at thaw.
  B: protocol change only; no formulation change.
  C: mechanical mixing of sealed vial; no contamination or formulation change.
     Risk: high-shear may affect shear-sensitive analytes (LDH, albumin) —
     documented in risk-assessment table and flagged in proposal text.
"""

import numpy as np
import csv
from pathlib import Path

from src.ripening_kinetics import (
    rate_constants as _rate_constants,
    R0_NM, K_LSW, VM, GAMMA, C_SAT_BULK,
    D_CA_22C, H_QUIESCENT, H_VORTEX, R_GAS,
)
from src.vial_simulation import (
    N_BATCHES, N_PER_BATCH, GLASS_REF, FILL_REF,
    NUCLEATION_MEDIAN_DAYS, NUCLEATION_SIGMA_LN_VIAL, NUCLEATION_SIGMA_LN_BATCH,
    BASE_K, THRESHOLD,
    _phase_fractions_analytical,
)

DATA_DIR    = Path(__file__).parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

N_VIALS       = N_BATCHES * N_PER_BATCH   # 10,000
THAW_TIME_MIN = 60.0
STORAGE_MONTHS = [6, 12]

# ── Scenario definitions ──────────────────────────────────────────────────────
#
# ph_storage:     pH in unfrozen cryo pool at −20°C
#   baseline 8.0  = 30% CO₂ outgassing (typical unsealed/loosely sealed vial)
#   +degas_50 7.9 = 50% residual CO₂ → less outgassing headroom
#   +degas_10 7.81= 10% residual CO₂ → ~0% outgassing (pKa₁ shift only)
#
# k_sigma:        σ of local_k_factor Gaussian
#   baseline 0.15 = ±15% spatial heterogeneity (uncontrolled freezer)
#   CRF_1C   0.08 = ±8% (1°C/min controlled rate, tighter ice morphology)
#   CRF_2C   0.05 = ±5% (2°C/min, smallest practical for amber glass)
#
# nuc_multiplier: scale factor on NUCLEATION_MEDIAN_DAYS
#   baseline  1.0 = nominal induction time
#   CRF adds: faster freeze → less time in intermediate-SI regime → longer delay
#   CRF_1C    1.8 (Pikal 2004: controlled freeze reduces nucleation frequency)
#   CRF_2C    2.5
#
# thaw_h:         Noyes-Whitney boundary layer thickness at thaw (m)
#   quiescent  10e-6 m  (baseline)
#   vortex_30s  2e-6 m  (literature: 1–5 µm in agitated dissolution)
#   vortex_60s  1.5e-6 m (more aggressive mixing, diminishing returns)

SCENARIOS: dict[str, dict] = {
    "baseline":    {"ph": 8.0,  "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": H_QUIESCENT},
    "+degas_50":   {"ph": 7.90, "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": H_QUIESCENT},
    "+degas_10":   {"ph": 7.81, "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": H_QUIESCENT},
    "+crf_1C":     {"ph": 8.0,  "k_sig": 0.08, "nuc_mult": 1.8, "thaw_h": H_QUIESCENT},
    "+crf_2C":     {"ph": 8.0,  "k_sig": 0.05, "nuc_mult": 2.5, "thaw_h": H_QUIESCENT},
    "+vortex_30s": {"ph": 8.0,  "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": H_VORTEX},
    "+vortex_60s": {"ph": 8.0,  "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": 1.5e-6},
    "+combined":   {"ph": 7.81, "k_sig": 0.05, "nuc_mult": 2.5, "thaw_h": H_VORTEX},
}


# ── Per-vial deficit calculator (scenario-aware) ──────────────────────────────

def _f_precip_scenario(albumin_g_dL: float, k_eff: float, ph: float) -> float:
    """Precipitated fraction at cryo state for given albumin, k, pH."""
    alb_pool = albumin_g_dL * k_eff
    alpha    = 1.0 / (1.0 + 0.25 * 10.0 ** (0.20 * (ph - 7.4)) * alb_pool)
    return min(1.0 - alpha, 0.98)


def compute_vial_deficit_scenario(
    storage_months: float,
    local_k_factor: float,
    nucleation_delay_days: float,
    particle_size_nm: float,
    storage_T_C: float,
    albumin_g_dL: float,
    cryo_purity: float,
    glass_density: float,
    fill_volume_mL: float,
    ph_storage: float,
    thaw_h: float,
    freezing_rate: float = 1.0,
) -> float:
    """Compute per-vial Ca deficit with scenario-specific pH and thaw protocol."""
    # Effective k
    fast_freeze_boost = 1.0 + 0.05 * np.log(max(freezing_rate, 0.1))
    k_eff = BASE_K * local_k_factor * cryo_purity * fast_freeze_boost

    # Nucleation delay
    surface_volume_ratio = fill_volume_mL / FILL_REF
    glass_factor         = (glass_density / GLASS_REF) * surface_volume_ratio
    effective_delay_days = nucleation_delay_days / max(glass_factor, 0.05)

    t_storage_days = storage_months * 30.4375
    t_eff_days     = max(0.0, t_storage_days - effective_delay_days)
    if t_eff_days < 1e-3:
        return 0.0

    # Rate constants (Arrhenius + pH-dependent)
    kk   = _rate_constants(storage_T_C, ph_storage)
    x_ACP, x_OCP, x_HAp = _phase_fractions_analytical(
        t_eff_days, kk["k_ACP_OCP"], kk["k_ACP_HAP"], kk["k_OCP_HAP"]
    )

    # Precipitated fraction (albumin at cryo state)
    fp = _f_precip_scenario(albumin_g_dL, k_eff, ph_storage)

    # Dissolution (Noyes-Whitney with scenario thaw_h)
    t_storage_h = t_storage_days * 24.0
    recovery    = 0.0

    for phase, frac, r0_override in [
        ("ACP", x_ACP, particle_size_nm),
        ("OCP", x_OCP, None),
        ("HAp", x_HAp, None),
    ]:
        if frac < 1e-6:
            continue
        r0   = r0_override if r0_override is not None else R0_NM[phase]
        r_nm = (r0**3 + K_LSW[phase] * t_storage_h) ** (1.0/3.0)
        r_m  = r_nm * 1e-9
        c_s  = C_SAT_BULK[phase] * np.exp(
            2.0 * GAMMA[phase] * VM[phase] / (r_m * R_GAS * 295.15)
        )
        c_prec   = 1.0 / (VM[phase] * 1000.0)
        h_eff    = max(r_m, thaw_h)
        lambda_p = min(D_CA_22C * (3.0 / r_m) / h_eff * (c_s / c_prec), 1.0)
        recovery += frac * (1.0 - np.exp(-lambda_p * THAW_TIME_MIN * 60.0))

    return fp * (1.0 - np.clip(recovery, 0.0, 1.0))


# ── Monte Carlo runner for a single scenario ──────────────────────────────────

def run_scenario(scenario_name: str, seed: int = 42) -> np.ndarray:
    """
    Return deficits array shape (N_VIALS, len(STORAGE_MONTHS)).
    Uses shared parameter sample (same seed) for all scenarios so that
    differences are due to interventions, not sampling noise.
    """
    sc = SCENARIOS[scenario_name]
    rng = np.random.default_rng(seed)

    # Sample batch parameters
    albumin       = rng.normal(5.5,  0.5,  N_BATCHES).clip(3.5,  7.5)
    cryo_purity   = rng.normal(1.0,  0.02, N_BATCHES).clip(0.94, 1.06)
    freezing_rate = np.exp(rng.normal(np.log(1.0), 0.3, N_BATCHES)).clip(0.3, 5.0)
    storage_T     = rng.normal(-20.0, 1.5, N_BATCHES).clip(-24.0, -16.0)
    nuc_batch_scale = np.exp(
        rng.normal(0.0, NUCLEATION_SIGMA_LN_BATCH, N_BATCHES)
    ).clip(0.05, 20.0)

    # Sample vial parameters
    batch_idx    = np.repeat(np.arange(N_BATCHES), N_PER_BATCH)
    local_k      = rng.normal(1.0, sc["k_sig"], N_VIALS).clip(0.5, 2.0)
    nuc_vial     = np.exp(
        rng.normal(
            np.log(NUCLEATION_MEDIAN_DAYS * sc["nuc_mult"]),
            NUCLEATION_SIGMA_LN_VIAL,
            N_VIALS,
        )
    )
    nuc_days     = nuc_vial * nuc_batch_scale[batch_idx]
    particle_nm  = np.exp(rng.normal(np.log(50.0), 0.5, N_VIALS)).clip(10.0, 300.0)
    glass_density = np.exp(rng.normal(np.log(0.5), 0.4, N_VIALS)).clip(0.05, 5.0)
    fill_volume  = rng.normal(7.5, 0.3, N_VIALS).clip(6.5, 8.5)

    deficits = np.zeros((N_VIALS, len(STORAGE_MONTHS)))

    for i in range(N_VIALS):
        bi = batch_idx[i]
        for j, sm in enumerate(STORAGE_MONTHS):
            deficits[i, j] = compute_vial_deficit_scenario(
                storage_months        = sm,
                local_k_factor        = local_k[i],
                nucleation_delay_days = nuc_days[i],
                particle_size_nm      = particle_nm[i],
                storage_T_C           = storage_T[bi],
                albumin_g_dL          = albumin[bi],
                cryo_purity           = cryo_purity[bi],
                glass_density         = glass_density[i],
                fill_volume_mL        = fill_volume[i],
                ph_storage            = sc["ph"],
                thaw_h                = sc["thaw_h"],
                freezing_rate         = freezing_rate[bi],
            )

    return deficits


# ── Compute all scenarios ─────────────────────────────────────────────────────

def run_all_scenarios(seed: int = 42) -> dict[str, np.ndarray]:
    results = {}
    for name in SCENARIOS:
        results[name] = run_scenario(name, seed=seed)
    return results


# ── Summary statistics ────────────────────────────────────────────────────────

def scenario_stats(all_results: dict[str, np.ndarray]) -> list[dict]:
    rows = []
    for name, deficits in all_results.items():
        row = {"scenario": name}
        for j, sm in enumerate(STORAGE_MONTHS):
            d = deficits[:, j]
            suffix = f"_{sm}mo"
            row[f"mean_deficit{suffix}_pct"]   = round(d.mean() * 100, 2)
            row[f"p95_deficit{suffix}_pct"]    = round(np.percentile(d, 95) * 100, 2)
            row[f"frac_above_4pct{suffix}"]    = round((d > THRESHOLD).mean(), 4)
        rows.append(row)
    return rows


def save_outcomes_csv(rows: list[dict], path: Path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {path}")


# ── Risk assessment table ─────────────────────────────────────────────────────

RISK_TABLE = [
    {
        "intervention":       "Baseline",
        "mechanism":          "No intervention",
        "pH_effect":          "pH 8.0 (30% CO₂ loss)",
        "IS_effect":          "None",
        "osmolality_effect":  "None",
        "metal_balance":      "Unaffected",
        "shear_risk":         "None",
        "complexity_1_5":     1,
        "capex_estimate":     "—",
        "ISO_13485":          "Compliant",
        "REACH_PFAS":         "Compliant",
        "validation_needed":  "—",
    },
    {
        "intervention":       "A: Pre-freeze degassing",
        "mechanism":          "Removes dissolved CO₂; prevents pH rise beyond pKa₁ shift",
        "pH_effect":          "pH 7.81 (reduces CO₂-driven rise by ~0.2 units)",
        "IS_effect":          "Negligible (CO₂ non-electrolyte)",
        "osmolality_effect":  "Negligible (<5 mOsm/kg)",
        "metal_balance":      "Unaffected (Ca, Mg, Na, K unchanged)",
        "shear_risk":         "None",
        "complexity_1_5":     2,
        "capex_estimate":     "Low (~$5-20k for inline vacuum degasser)",
        "ISO_13485":          "Compliant (gas removal only; chemical composition unchanged at thaw)",
        "REACH_PFAS":         "Compliant",
        "validation_needed":  ("pH measurement degassed vs control batches; "
                               "residual CO₂ by GC-headspace; "
                               "Ca recovery comparison 6-month storage"),
    },
    {
        "intervention":       "B: Controlled-rate freezing (1°C/min)",
        "mechanism":          ("Reduces spatial k heterogeneity; smaller ice crystals; "
                               "less time at intermediate SI → longer nucleation delay"),
        "pH_effect":          "None (pH trajectory unchanged)",
        "IS_effect":          "None",
        "osmolality_effect":  "None",
        "metal_balance":      "Unaffected",
        "shear_risk":         "None",
        "complexity_1_5":     4,
        "capex_estimate":     "High (~$30-100k per controlled-rate freezer unit)",
        "ISO_13485":          "Compliant (validated process change under ICH Q14)",
        "REACH_PFAS":         "Compliant",
        "validation_needed":  ("Cryo-SEM ice morphology comparison; "
                               "thermal trajectory monitoring (IQ/OQ/PQ); "
                               "6-month comparative storage study"),
    },
    {
        "intervention":       "C: Modified thaw (30-s vortex at 30 min)",
        "mechanism":          ("Reduces Noyes-Whitney boundary layer ~5×; "
                               "HAp microcrystals dissolve in minutes instead of hours"),
        "pH_effect":          "None at 60 min; negligible drift (<0.05 pH) from CO₂ equilibration",
        "IS_effect":          "None",
        "osmolality_effect":  "None",
        "metal_balance":      "Unaffected for Ca/Mg; RISK: LDH activity may decrease "
                               "≤15% at 1500 rpm 30-s; albumin aggregation risk low but non-zero",
        "shear_risk":         ("MODERATE — shear-sensitive analytes: LDH, CK, "
                               "protein-bound hormones; MUST validate empirically"),
        "complexity_1_5":     2,
        "capex_estimate":     "Minimal (<$500; vortex mixer already present in most labs)",
        "ISO_13485":          "Compliant; requires SOP change + analytical validation per ICH Q2(R2)",
        "REACH_PFAS":         "Compliant",
        "validation_needed":  ("Empirical dose-response: vortex duration vs Ca recovery; "
                               "analyte stability panel (LDH, CK, albumin, cortisol, T4) "
                               "after vortex; no-spin control comparison"),
    },
    {
        "intervention":       "Combined (A+B+C)",
        "mechanism":          "All three interventions applied together",
        "pH_effect":          "pH 7.81 (degassing controls pH)",
        "IS_effect":          "None",
        "osmolality_effect":  "Negligible",
        "metal_balance":      "Unaffected for Ca/Mg; shear risk as per C",
        "shear_risk":         "Moderate (same as C)",
        "complexity_1_5":     5,
        "capex_estimate":     "High (dominated by controlled-rate freezer)",
        "ISO_13485":          "Compliant; full validation required",
        "REACH_PFAS":         "Compliant",
        "validation_needed":  ("Combination of individual validation studies above; "
                               "full design-of-experiments recommended"),
    },
]


def save_risk_table(path: Path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(RISK_TABLE[0].keys()))
        w.writeheader()
        w.writerows(RISK_TABLE)
    print(f"  Saved: {path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main_interventions():
    print("Module 7: Intervention modeling")
    print("=" * 60)
    print(f"Scenarios: {list(SCENARIOS.keys())}")
    print(f"Storage evaluated at: {STORAGE_MONTHS} months")
    print()

    print("Running all scenarios …")
    all_results = run_all_scenarios(seed=42)

    rows = scenario_stats(all_results)

    # Print summary table
    print(f"{'Scenario':<18}  {'6mo frac>4%':>12}  {'12mo frac>4%':>13}  "
          f"{'12mo mean%':>10}  {'12mo p95%':>9}")
    for row in rows:
        print(f"  {row['scenario']:<16}  {row['frac_above_4pct_6mo']:>12.3f}  "
              f"{row['frac_above_4pct_12mo']:>13.3f}  "
              f"{row['mean_deficit_12mo_pct']:>10.1f}  "
              f"{row['p95_deficit_12mo_pct']:>9.1f}")

    save_outcomes_csv(rows, DATA_DIR / "module7_intervention_outcomes.csv")
    save_risk_table(DATA_DIR / "module7_risk_assessment.csv")

    return all_results, rows


if __name__ == "__main__":
    main_interventions()
