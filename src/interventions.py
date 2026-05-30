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
    D_CA_22C, D_CA_4C, H_QUIESCENT, H_VORTEX, R_GAS, F_PRECIP, ACP_AGGREGATE_NM,
    nucleation_temp_factor,
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
# PRIMARY BASELINE is the SEALED vial (pH 7.81, no CO₂ outgassing) — the most
# defensible physical state for a closed QC vial. `loose_seal` (pH 8.0) is kept
# as an explicit RISK scenario: if the closure leaks and CO₂ escapes, pH rises
# and the problem gets worse — which is itself the rationale for degassing /
# tight sealing. (Earlier versions used pH 8.0 as the baseline in this module
# while Module 6 used 7.81; that mismatch produced the 42% vs 58% discrepancy.
# Both are now reported as named, distinct states.)
#
# ph_storage:     pH in unfrozen cryo pool at −20°C
#   baseline   7.81 = sealed vial (pKa shift + ionic strength only)
#   loose_seal 8.0  = CO₂ outgassing through a leaking closure (risk state)
#
# k_sigma:        σ of local_k_factor Gaussian
#   baseline 0.15 = ±15% spatial heterogeneity (uncontrolled freezer)
#   CRF_2C   0.05 = ±5% (2°C/min controlled-rate freezing)
#
# nuc_multiplier: scale factor on NUCLEATION_MEDIAN_DAYS
#   baseline 1.0; CRF 2.5 (faster freeze → longer nucleation delay; Pikal 2004)
#
# thaw_h:         Noyes-Whitney boundary layer thickness at thaw (m)
#   quiescent  10e-6 m   (baseline, no mixing)
#   vortex_30s  2e-6 m   (5× mass transfer; literature 1–5 µm agitated)
#   vortex_60s  1.5e-6 m (≈7× mass transfer; diminishing returns)
#   combined_plus 1.0e-6 m (double-pulse vortex → 10× mass transfer)
#
# thaw_min:       Measurement window (min): 60 standard; 90 combined_plus
#                 (60 RT + 30 cold soak); 2880 = extended-mixing reference (48 h).
# d_factor:       Diffusivity multiplier (1.0 = 22°C; ≈0.60 = 4°C, Stokes-Einstein)

_D_FACTOR_4C = D_CA_4C / D_CA_22C   # ≈ 0.60 at 4°C vs 22°C

SCENARIOS: dict[str, dict] = {
    "baseline":    {"ph": 7.81, "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": H_QUIESCENT},
    # Risk state (not an intervention): leaking closure → CO₂ loss → higher pH.
    "loose_seal":  {"ph": 8.0,  "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": H_QUIESCENT},
    "+vortex_30s": {"ph": 7.81, "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": H_VORTEX},
    "+vortex_60s": {"ph": 7.81, "k_sig": 0.15, "nuc_mult": 1.0, "thaw_h": 1.5e-6},
    "+crf_2C":     {"ph": 7.81, "k_sig": 0.05, "nuc_mult": 2.5, "thaw_h": H_QUIESCENT},
    "+combined":   {"ph": 7.81, "k_sig": 0.05, "nuc_mult": 2.5, "thaw_h": H_VORTEX},
    # Double-pulse vortex: 30 s at 5 min + 60 s at 25 min → h=1 µm; plus 30-min
    # cold soak at 2-8°C before measurement (total measurement window 90 min).
    "+combined_plus": {
        "ph": 7.81, "k_sig": 0.05, "nuc_mult": 2.5,
        "thaw_h": 1.0e-6, "thaw_min": 90.0,
    },
    # Extended-mixing / cold-equilibration reference: long quiescent soak at
    # 2-8°C. Confirms that, given enough time, the deficit fully reverses —
    # consistent with the Seeker's report that mixing/standing resolves it.
    "+extended_mixing": {
        "ph": 7.81, "k_sig": 0.15, "nuc_mult": 1.0,
        "thaw_h": H_QUIESCENT, "thaw_min": 2880.0, "d_factor": _D_FACTOR_4C,
    },
    # PREVENTION (root cause): deep-frozen / vitrified storage at ≤ −80°C. Below
    # the freeze-concentrate glass transition the pool is immobile, nucleation is
    # arrested, and the precipitate never forms — so even a standard quiescent
    # thaw reads correct. No mixing needed; the fix is upstream, at storage.
    "+deep_freeze": {
        "ph": 7.81, "k_sig": 0.15, "nuc_mult": 1.0,
        "thaw_h": H_QUIESCENT, "storage_T_mean": -80.0,
    },
}

# Scenarios that are genuine interventions (exclude baseline + the risk state).
INTERVENTION_NAMES = [n for n in SCENARIOS if n not in ("baseline", "loose_seal")]


# ── Per-vial deficit calculator (scenario-aware) ──────────────────────────────

def _f_precip_scenario(albumin_g_dL: float = 0.0, k_eff: float = 0.0,
                       ph: float = 7.81) -> float:
    """Precipitated fraction: mass-balance constant (see F_PRECIP note).

    Arguments retained for compatibility but unused — the precipitated fraction
    is not an albumin-binding quantity.
    """
    return F_PRECIP


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
    thaw_min: float = 60.0,
    d_factor: float = 1.0,
    freezing_rate: float = 1.0,
) -> float:
    """Compute per-vial Ca deficit with scenario-specific pH and thaw protocol.

    thaw_min: measurement window in minutes (60 = standard; 2880 = 48-h soak).
    d_factor: diffusivity multiplier for thaw temperature (1.0 = 22°C, ~0.60 = 4°C).
    """
    # Effective k
    fast_freeze_boost = 1.0 + 0.05 * np.log(max(freezing_rate, 0.1))
    k_eff = BASE_K * local_k_factor * cryo_purity * fast_freeze_boost

    # Nucleation delay
    surface_volume_ratio = fill_volume_mL / FILL_REF
    glass_factor         = (glass_density / GLASS_REF) * surface_volume_ratio
    effective_delay_days = nucleation_delay_days / max(glass_factor, 0.05)

    # Supersaturation control of nucleation: a higher pool pH means more
    # supersaturation and a shorter induction time (classical nucleation theory,
    # monotonic heuristic). This is how degassing/tight sealing (which keep pH
    # low) and CO₂ outgassing (which raises pH) act — on the AFFECTED FRACTION
    # via nucleation, not on the per-vial deficit magnitude.
    ph_super_factor       = 10.0 ** (2.0 * (ph_storage - 7.81))
    effective_delay_days /= ph_super_factor

    # PREVENTION lever: colder storage → much more viscous pool → longer
    # induction. Below Tg' (~−50°C) the pool vitrifies and nucleation is
    # arrested → deep-frozen storage prevents the precipitate from forming.
    effective_delay_days *= nucleation_temp_factor(storage_T_C)

    t_storage_days = storage_months * 30.4375
    t_eff_days     = max(0.0, t_storage_days - effective_delay_days)
    if t_eff_days < 1e-3:
        return 0.0

    # Rate constants (kept for reference; ripening to HAp is suppressed at the
    # corrected pool viscosity, so the precipitate stays essentially all ACP).
    kk   = _rate_constants(storage_T_C, ph_storage)
    x_ACP, x_OCP, x_HAp = _phase_fractions_analytical(
        t_eff_days, kk["k_ACP_OCP"], kk["k_ACP_HAP"], kk["k_OCP_HAP"]
    )

    # Precipitated fraction (mass-balance constant)
    fp = _f_precip_scenario()

    # Re-dispersion of the wall-bound amorphous precipitate. The thaw protocol
    # sets the diffusion boundary layer thaw_h (mixing → thin → fast), the window
    # thaw_min, and the diffusivity (cold soak). Deficit = precipitated fraction
    # × fraction NOT redispersed.
    D_thaw = D_CA_22C * d_factor
    r_m    = ACP_AGGREGATE_NM * 1e-9
    c_s    = C_SAT_BULK["ACP"] * np.exp(2.0*GAMMA["ACP"]*VM["ACP"]/(r_m*R_GAS*295.15))
    c_prec = 1.0 / (VM["ACP"] * 1000.0)
    h_eff  = max(r_m, thaw_h)
    lam    = min(D_thaw * (3.0 / r_m) / h_eff * (c_s / c_prec), 1.0)
    recovery = 1.0 - np.exp(-lam * thaw_min * 60.0)

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
    T_mean        = sc.get("storage_T_mean", -20.0)
    storage_T     = rng.normal(T_mean, 1.5, N_BATCHES).clip(T_mean - 4.0, T_mean + 4.0)
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
                thaw_min              = sc.get("thaw_min", 60.0),
                d_factor              = sc.get("d_factor", 1.0),
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
            row[f"frac_with_deficit{suffix}"]  = round((d > THRESHOLD).mean(), 4)
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
    {
        "intervention":       "D: Combined+ (A+B+double-pulse vortex+cold soak)",
        "mechanism":          ("As Combined, but thaw protocol: 30-s vortex at 5 min, "
                               "60-s vortex at 25 min, then 30-min equilibration at 2-8°C "
                               "before measurement. Two vortex pulses reduce boundary layer "
                               "to ~1 µm; cold soak adds dissolution time."),
        "pH_effect":          "pH 7.81 (degassing controls pH)",
        "IS_effect":          "None",
        "osmolality_effect":  "Negligible",
        "metal_balance":      "Unaffected for Ca/Mg; double-vortex shear risk slightly higher",
        "shear_risk":         ("MODERATE-HIGH — two vortex pulses at 1500 rpm; "
                               "MUST validate analyte stability panel"),
        "complexity_1_5":     5,
        "capex_estimate":     "High (as Combined; vortex already present)",
        "ISO_13485":          "Compliant; requires updated SOP for thaw protocol",
        "REACH_PFAS":         "Compliant",
        "validation_needed":  ("All Combined validation studies; dose-response for "
                               "double-pulse vortex; analyte stability at 2-8°C cold soak"),
    },
    {
        "intervention":       "Seeker workaround: 48-h cold equilibration at 2-8°C",
        "mechanism":          ("Extended equilibration at 2-8°C provides sufficient time "
                               "for HAp microcrystals to dissolve toward solubility limit. "
                               "At 4°C, D(Ca²⁺) ≈ 60% of 22°C value, but 48 h = 2880 min "
                               "gives ~28× more dissolution time than standard 60-min thaw."),
        "pH_effect":          "Unchanged (pH 8.0 baseline)",
        "IS_effect":          "None",
        "osmolality_effect":  "None",
        "metal_balance":      "Unaffected",
        "shear_risk":         "None (no vortex)",
        "complexity_1_5":     1,
        "capex_estimate":     "Zero (cold room already present); only SOP change",
        "ISO_13485":          "Compliant; SOP change only; no new equipment",
        "REACH_PFAS":         "Compliant",
        "validation_needed":  ("Turnaround time 48 h is operationally impractical for "
                               "routine use; 24 h may suffice; validate minimum soak time; "
                               "confirm stable Ca at 48 h vs 6 h cold soak"),
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
    print(f"{'Scenario':<18}  {'6mo frac w/def':>14}  {'12mo frac w/def':>15}  "
          f"{'12mo mean%':>10}")
    for row in rows:
        print(f"  {row['scenario']:<16}  {row['frac_with_deficit_6mo']:>14.3f}  "
              f"{row['frac_with_deficit_12mo']:>15.3f}  "
              f"{row['mean_deficit_12mo_pct']:>10.1f}")

    save_outcomes_csv(rows, DATA_DIR / "module7_intervention_outcomes.csv")
    save_risk_table(DATA_DIR / "module7_risk_assessment.csv")

    return all_results, rows


if __name__ == "__main__":
    main_interventions()
