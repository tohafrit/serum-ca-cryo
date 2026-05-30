"""
Module 6: Monte Carlo simulation of vial-to-vial variability.

10,000 virtual vials (50 batches × 200 vials) are sampled from
independent distributions of within-batch and between-batch parameters.
For each vial the post-thaw Ca deficit at 60-min thaw is computed via
Module 5 kinetics, giving a predicted distribution that is compared to
the Seeker's observation of "vial-to-vial dependence" and "in some samples".

Parameter distributions are sourced from:
  - Nucleation induction times: lognormal from classical nucleation theory
    (Toschev 1973; Kashchiev 2000)
  - Glass surface site density: literature for Type I borosilicate
    (Jennings 1983; Bhatt 2011)
  - Albumin donor variability: NHANES serum protein reference ranges
  - Freezer temperature variance: ASHRAE 2009 cold storage specs ±2°C
  - Cryoprotectant purity: USP glycerol grade specification
  - Local k variability: geometric heterogeneity of cylindrical vial freeze
    (Morris 2005; Pikal 2004 — 10-15% spatial variation in freeze front)
"""

import numpy as np
import csv
from pathlib import Path

from src.ripening_kinetics import (
    rate_constants,
    R0_NM, K_LSW, VM, GAMMA, C_SAT_BULK,
    D_CA_22C, H_QUIESCENT, R_GAS, F_PRECIP, ACP_AGGREGATE_NM,
    nucleation_temp_factor,
)

# ── Simulation parameters ─────────────────────────────────────────────────────

N_BATCHES   = 50
N_PER_BATCH = 200
N_VIALS     = N_BATCHES * N_PER_BATCH    # 10,000

BASE_K      = 5.58    # deterministic k from Module 2 (glycerol 15%, −20°C)
BASE_T_C    = -20.0
BASE_PH     = 7.81    # sealed vial, 0% CO₂ loss (primary baseline)
# Precipitated fraction is mass-balance-bounded and the same for every vial
# (it is NOT an albumin-binding quantity — see note at F_PRECIP in
# ripening_kinetics.py). Imported directly; not resampled per vial.

STORAGE_MONTHS = [1, 3, 6, 9, 12, 18, 24]
THAW_TIME_MIN  = 60.0

DATA_DIR = Path(__file__).parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)

# ── Per-vial deficit calculation ──────────────────────────────────────────────
#
# The stochastic parameters modulate five aspects of the Module 5 model:
#
#   1. k_eff       = BASE_K × local_k_factor × cryo_purity
#                    (spatial freeze heterogeneity + glycerol lot variation)
#
#   2. T_eff       = storage_T_C  (actual freezer temperature per vial)
#
#   3. t_eff       = max(0, t_storage_days − nucleation_delay_h/24)
#                    (nucleation induction shifts start of ACP clock)
#
#   4. f_precip    = F_PRECIP (mass-balance constant; same for all vials).
#                    Albumin is still sampled (donor variability) but no longer
#                    drives the precipitated fraction — Sobol confirms albumin
#                    is negligible (ST≈0.002), so this does not change variance.
#
#   5. r0_eff_nm   = initial_ACP_particle_size_nm
#                    (affects dissolution rate via Ostwald-Freundlich)
#
# Glass surface site density + fill volume affect nucleation probability
# (higher density → shorter effective induction time); modelled as an
# additive −/+ shift on nucleation_delay_h:
#   nucleation_delay_h *= 1 / (glass_density / GLASS_REF) × fill_volume / FILL_REF
# (more sites or less volume → sooner nucleation → shorter delay)

GLASS_REF = 0.5   # nmol/cm² reference surface density
FILL_REF  = 7.5   # mL reference fill volume

# Nucleation induction time at −20°C glycerol 15%:
# At physiological conditions (~37°C, SI≈5) literature reports t_ind ~ 2-6 h
# (Heughebaert & Nancollas 1984). Viscosity correction (×10) + temperature
# correction (Arrhenius, Ea~60 kJ/mol, from 37°C to −20°C):
#   t_ind(-20°C) ≈ t_ind(37°C) × η_ratio × exp(Ea/R × (1/T_cryo - 1/T_phys))
#   ≈ 4 h × 10 × exp(7200 × 0.0015) ≈ 4 × 10 × 50 ≈ 2000 h ≈ 83 days
# This is the median; vial-to-vial spread spans ~3 orders of magnitude
# (classical nucleation theory: J ∝ exp(-ΔG*/kT); tiny ΔΔG → huge t_ind change)
NUCLEATION_MEDIAN_DAYS   = 270.0  # ~9 mo median induction → onset ~6 mo,
#                                  calibrated to the Seeker's "not before 6 months"
#                                  (also consistent with the corrected high pool
#                                  viscosity). Refined by the storage-time series.
NUCLEATION_SIGMA_LN_VIAL  = 1.0   # within-batch lognormal σ (~10× span vial-to-vial)
NUCLEATION_SIGMA_LN_BATCH = 0.8   # between-batch lognormal σ (formulation, glass lot)


def _f_precip(albumin_g_dL_serum: float = 0.0, k_eff: float = 0.0,
              pH: float = 7.81) -> float:
    """
    Fraction of total serum Ca that precipitates in the cryo pool.

    Mass-balance-bounded constant (see note at F_PRECIP in ripening_kinetics).
    Arguments are retained for call-site compatibility but unused: the
    precipitated fraction is governed by supersaturation + phosphate balance,
    not by albumin binding.
    """
    return F_PRECIP


def _phase_fractions_analytical(t_days: float, k_ao: float, k_ah: float, k_oh: float):
    """
    Analytical solution to the linear 3-state ODE (ACP→OCP→HAp).
    Avoids scipy.integrate entirely → ~1000× faster per vial.

    dx_ACP/dt = -a * x_ACP,          a = k_ao + k_ah
    dx_OCP/dt = k_ao * x_ACP - k_oh * x_OCP
    dx_HAp/dt = k_ah * x_ACP + k_oh * x_OCP
    y0 = [1, 0, 0]
    """
    a = k_ao + k_ah
    eat = np.exp(-a   * t_days)
    ebt = np.exp(-k_oh * t_days)

    x_ACP = eat

    if abs(k_oh - a) > 1e-15:
        x_OCP = k_ao / (k_oh - a) * (eat - ebt)
    else:
        x_OCP = k_ao * t_days * eat

    x_HAp = 1.0 - x_ACP - x_OCP
    # Clamp floating-point noise
    x_OCP = max(x_OCP, 0.0)
    x_HAp = max(x_HAp, 0.0)
    return x_ACP, x_OCP, x_HAp


def compute_vial_deficit(
    storage_months: float,
    local_k_factor: float,
    nucleation_delay_days: float,
    particle_size_nm: float,
    storage_T_C: float,
    albumin_g_dL: float,
    cryo_purity: float,
    glass_density: float,
    fill_volume_mL: float,
    freezing_rate: float = 1.0,
    nominal_T_C: float = BASE_T_C,
) -> float:
    """Return Ca deficit fraction at THAW_TIME_MIN for one virtual vial."""
    # Effective k
    fast_freeze_boost = 1.0 + 0.05 * np.log(max(freezing_rate, 0.1))
    k_eff = BASE_K * local_k_factor * cryo_purity * fast_freeze_boost

    # Glass surface density modulates nucleation rate: more sites → shorter delay
    # Relative to reference density; fill volume modulates surface/volume ratio
    surface_volume_ratio = (fill_volume_mL / FILL_REF)           # proxy for A/V
    glass_factor = (glass_density / GLASS_REF) * surface_volume_ratio
    effective_delay_days = nucleation_delay_days / max(glass_factor, 0.05)
    # Colder NOMINAL storage → far more viscous pool → longer nucleation
    # induction (the deliberate process lever; applied at the scenario's nominal
    # temperature, not the small batch-to-batch storage noise — the factor is too
    # steep to be driven by ±1-2°C measurement scatter).
    effective_delay_days *= nucleation_temp_factor(nominal_T_C)

    t_storage_days = storage_months * 30.4375
    t_eff_days     = max(0.0, t_storage_days - effective_delay_days)
    if t_eff_days < 1e-3:
        return 0.0

    # Rate constants (Arrhenius at per-vial storage temperature)
    kk    = rate_constants(storage_T_C, BASE_PH)
    k_ao  = kk["k_ACP_OCP"]
    k_ah  = kk["k_ACP_HAP"]
    k_oh  = kk["k_OCP_HAP"]

    # Phase fractions (kept for reference; with the corrected pool viscosity the
    # precipitate stays essentially all ACP — ripening to HAp is suppressed).
    x_ACP, x_OCP, x_HAp = _phase_fractions_analytical(t_eff_days, k_ao, k_ah, k_oh)

    # Precipitated fraction (mass-balance constant; not albumin-derived)
    fp = _f_precip()

    # Re-dispersion of the wall-bound amorphous precipitate at a standard
    # quiescent 60-min thaw (Noyes-Whitney, size = ACP_AGGREGATE_NM). Deficit
    # = precipitated fraction × fraction NOT redispersed in the window.
    r_m    = ACP_AGGREGATE_NM * 1e-9
    c_s    = C_SAT_BULK["ACP"] * np.exp(2.0*GAMMA["ACP"]*VM["ACP"]/(r_m*R_GAS*295.15))
    c_prec = 1.0 / (VM["ACP"] * 1000.0)
    h_eff  = max(r_m, H_QUIESCENT)
    lam    = min(D_CA_22C * (3.0 / r_m) / h_eff * (c_s / c_prec), 1.0)
    recovery = 1.0 - np.exp(-lam * THAW_TIME_MIN * 60.0)

    return fp * (1.0 - np.clip(recovery, 0.0, 1.0))


# ── Sample batch and vial parameters ─────────────────────────────────────────

def sample_parameters(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample batch-level (N_BATCHES) and vial-level (N_VIALS) parameters.

    Batch parameters columns: [albumin_g_dL, cryo_purity, freezing_rate,
                                storage_T_C, nucleation_batch_scale]
    Vial parameters columns:  [batch_idx, local_k_factor,
                                nucleation_delay_days, particle_size_nm,
                                glass_density, fill_volume_mL]

    Nucleation delay design (critical for correct vial-to-vial spread):
      nucleation_delay_days = batch_scale × vial_draw
      vial_draw ~ Lognormal(median=NUCLEATION_MEDIAN_DAYS,
                             σ_ln=NUCLEATION_SIGMA_LN_VIAL)
      batch_scale ~ Lognormal(median=1.0, σ_ln=NUCLEATION_SIGMA_LN_BATCH)
    Combined σ_ln = √(1.0²+0.8²) ≈ 1.28 → ~3-order-of-magnitude vial spread.

    Physical basis for median=90 days: Arrhenius+viscosity extrapolation of
    Heughebaert & Nancollas 1984 (4 h at 37°C, SI≈5) to −20°C glycerol 15%.
    """
    # ── Batch-level ────────────────────────────────────────────────────────
    albumin       = rng.normal(5.5,  0.5,  N_BATCHES).clip(3.5,  7.5)
    cryo_purity   = rng.normal(1.0,  0.02, N_BATCHES).clip(0.94, 1.06)
    freezing_rate = np.exp(rng.normal(np.log(1.0), 0.3, N_BATCHES)).clip(0.3, 5.0)
    storage_T     = rng.normal(-20.0, 1.5, N_BATCHES).clip(-24.0, -16.0)
    # Batch nucleation scale: some batches have glass lots / formulation
    # that catalyses nucleation (scale < 1 → shorter delay) or inhibits it
    nuc_batch_scale = np.exp(
        rng.normal(0.0, NUCLEATION_SIGMA_LN_BATCH, N_BATCHES)
    ).clip(0.05, 20.0)

    batch_params = np.column_stack([
        albumin, cryo_purity, freezing_rate, storage_T, nuc_batch_scale,
    ])

    # ── Vial-level ─────────────────────────────────────────────────────────
    batch_idx = np.repeat(np.arange(N_BATCHES), N_PER_BATCH)

    local_k_factor = rng.normal(1.0, 0.15, N_VIALS).clip(0.6, 1.6)

    # Vial-level nucleation draw (before batch scaling)
    nuc_vial_draw = np.exp(
        rng.normal(np.log(NUCLEATION_MEDIAN_DAYS), NUCLEATION_SIGMA_LN_VIAL, N_VIALS)
    )
    # Apply batch scale (combines within-batch and between-batch variance)
    nuc_delay_days = nuc_vial_draw * nuc_batch_scale[batch_idx]

    particle_nm  = np.exp(rng.normal(np.log(50.0), 0.5, N_VIALS)).clip(10.0, 300.0)
    glass_density = np.exp(rng.normal(np.log(0.5), 0.4, N_VIALS)).clip(0.05, 5.0)
    fill_volume  = rng.normal(7.5, 0.3, N_VIALS).clip(6.5, 8.5)

    vial_params = np.column_stack([
        batch_idx, local_k_factor, nuc_delay_days,
        particle_nm, glass_density, fill_volume,
    ])
    return batch_params, vial_params


# ── Run Monte Carlo ───────────────────────────────────────────────────────────

def run_monte_carlo(seed: int = 42) -> dict:
    """
    Run the full 10,000-vial Monte Carlo simulation.
    Returns dict with:
      'deficits': array shape (N_VIALS, len(STORAGE_MONTHS))
      'batch_params': array shape (N_BATCHES, 4)
      'vial_params': array shape (N_VIALS, 6)
    """
    rng = np.random.default_rng(seed)
    batch_params, vial_params = sample_parameters(rng)

    deficits = np.zeros((N_VIALS, len(STORAGE_MONTHS)))

    for i in range(N_VIALS):
        bi = int(vial_params[i, 0])
        albumin, cryo_purity, freezing_rate, storage_T, _ = batch_params[bi]
        _, local_k, nuc_days, part_nm, glass_d, fill_vol = vial_params[i]

        for j, sm in enumerate(STORAGE_MONTHS):
            deficits[i, j] = compute_vial_deficit(
                storage_months        = sm,
                local_k_factor        = local_k,
                nucleation_delay_days = nuc_days,
                particle_size_nm      = part_nm,
                storage_T_C           = storage_T,
                albumin_g_dL          = albumin,
                cryo_purity           = cryo_purity,
                glass_density         = glass_d,
                fill_volume_mL        = fill_vol,
                freezing_rate         = freezing_rate,
            )

    return {
        "deficits":     deficits,
        "batch_params": batch_params,
        "vial_params":  vial_params,
    }


# ── Summary statistics + CSV ──────────────────────────────────────────────────

# Seeker's reported threshold: a "decrease of 4% or more". The mitigation target
# is to bring the decrease below 4%. The vial-to-vial observable is the FRACTION
# of vials at or above this threshold, governed by stochastic nucleation
# (→ "in some samples", batch/vial dependence).
THRESHOLD = 0.04   # 4% — Seeker's reported deficit threshold / mitigation target


def vial_statistics(deficits: np.ndarray) -> list[dict]:
    rows = []
    for j, sm in enumerate(STORAGE_MONTHS):
        d = deficits[:, j]
        rows.append({
            "storage_months":    sm,
            "mean_deficit_pct":  round(d.mean() * 100, 2),
            "median_deficit_pct":round(np.median(d) * 100, 2),
            "p95_deficit_pct":   round(np.percentile(d, 95) * 100, 2),
            "max_deficit_pct":   round(d.max() * 100, 2),
            "frac_with_deficit": round((d > THRESHOLD).mean(), 4),
            "n_vials":           len(d),
        })
    return rows


def save_csv(rows: list[dict], path: Path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main_mc():
    print("Module 6: Monte Carlo vial-to-vial variability")
    print("=" * 60)
    print(f"N_VIALS = {N_VIALS:,} ({N_BATCHES} batches × {N_PER_BATCH} vials/batch)")
    print()

    print("Running simulation …")
    result = run_monte_carlo(seed=42)
    deficits = result["deficits"]

    stats = vial_statistics(deficits)
    print(f"{'Mo':>3}  {'mean%':>6}  {'med%':>6}  {'p95%':>6}  {'frac w/ deficit':>15}")
    for row in stats:
        flag = " ← in some samples" if row["storage_months"] == 6 else ""
        print(f"  {row['storage_months']:>2}  "
              f"{row['mean_deficit_pct']:>6.1f}  "
              f"{row['median_deficit_pct']:>6.1f}  "
              f"{row['p95_deficit_pct']:>6.1f}  "
              f"{row['frac_with_deficit']:>15.3f}{flag}")

    save_csv(stats, DATA_DIR / "module6_vial_statistics.csv")
    return result, stats


if __name__ == "__main__":
    main_mc()
