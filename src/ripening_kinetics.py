"""
Module 5: Ostwald ripening kinetics and post-thaw redissolution.

Part A — kinetics during storage at −20°C:
  Three-state ODE: ACP → OCP → HAp
  Rate constants from Arrhenius extrapolation of literature data.

Part B — post-thaw redissolution:
  Noyes-Whitney with Ostwald-Freundlich size correction.
  Reports Ca-recovery vs time for four protocols × four storage durations.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import norm
from pathlib import Path
import csv

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR    = Path(__file__).parent.parent / "data"
FIGURES_DIR = Path(__file__).parent.parent / "figures"
DATA_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# PART A: ACP → OCP → HAp transformation kinetics at −20°C
# ─────────────────────────────────────────────────────────────────────────────
#
# Mechanism: solution-mediated dissolution-reprecipitation (not solid-state).
# ACP dissolves into the (supersaturated) unfrozen pool, OCP/HAp nucleate and
# grow.  At −20°C the unfrozen pool is ~18 wt% of the total water volume
# (glycerol 15%, k≈5.58).  Everything takes place in that thin film.
#
# Rate constants calibrated to:
#   Boskey & Posner 1973 — ACP→HAp t½ ≈ 24 h at 25°C, neutral pH
#   Eanes & Posner 1970  — ACP→HAp t½ ≈ 6 h at 37°C, pH 7.4
#   Combes & Rey 2010    — ACP kinetically stable for months at <0°C
#   Christoffersen 1989  — Ea (ACP→OCP) ≈ 65 kJ/mol
#   Dorozhkin 2010       — OCP→HAp slower than ACP→OCP; Ea ≈ 80 kJ/mol
#
# Arrhenius extrapolation:
#   k(T) = k(T_ref) × exp(-Ea/R × (1/T − 1/T_ref))
#
# pH scaling:
#   ACP→OCP rate increases ~10× per pH unit above 7 (Christoffersen 1989);
#   primary scenario uses pH=7.81 (0% CO₂ loss at −20°C).

R_GAS   = 8.314   # J/(mol·K)
T_REF_K = 298.15  # 25°C reference

# Literature rate constants at 25°C, neutral pH (first-order, per day)
# ACP→OCP: inferred so that k_ACP→OCP + k_ACP→HAp gives t½(ACP)≈24 h at 25°C
# OCP→HAp: slower; Dorozhkin 2010 gives t½(OCP) weeks at 25°C → k≈0.05 /day
K_ACP_OCP_25C = 0.60   # /day at 25°C, pH 7
K_ACP_HAP_25C = 0.10   # /day at 25°C, pH 7 (minor direct path)
K_OCP_HAP_25C = 0.05   # /day at 25°C, pH 7

EA_ACP_OCP = 65_000.0  # J/mol (Christoffersen 1989)
EA_ACP_HAP = 65_000.0  # J/mol (same mechanism, assume equal)
EA_OCP_HAP = 80_000.0  # J/mol (Dorozhkin 2010; solid-state slower)

# Glycerol-15% viscosity at −20°C ≈ 10× water at 0°C (Sheely 1932).
# Solution-mediated ripening is diffusion-limited → apply Stokes-Einstein:
#   k_diff ∝ D ∝ T/η  →  scale = (T_cryo/η_cryo) / (T_ref/η_ref)
# η_water(25°C) = 0.89 mPa·s; η_glycerol15%(−20°C) ≈ 9 mPa·s
VISC_SCALE = (253.15 / 9.0) / (298.15 / 0.89)   # ≈ 0.094

# pH effect on ACP→OCP: 10× per pH unit above 7.0 (sigmoidal, approximated)
def _ph_factor(pH: float) -> float:
    """Rate multiplier relative to pH 7.0 (Christoffersen 1989 scaling)."""
    return 10.0 ** (pH - 7.0)


def rate_constants(T_celsius: float, pH: float) -> dict:
    """
    Return first-order rate constants (units: 1/day) for transformation
    at given temperature and pH, with Arrhenius and viscosity corrections.
    """
    T_K    = T_celsius + 273.15
    arr    = lambda Ea: np.exp(-Ea / R_GAS * (1.0/T_K - 1.0/T_REF_K))
    ph_fac = _ph_factor(pH)

    k_ao = K_ACP_OCP_25C * arr(EA_ACP_OCP) * VISC_SCALE * ph_fac
    k_ah = K_ACP_HAP_25C * arr(EA_ACP_HAP) * VISC_SCALE * ph_fac
    k_oh = K_OCP_HAP_25C * arr(EA_OCP_HAP) * VISC_SCALE

    return {"k_ACP_OCP": k_ao, "k_ACP_HAP": k_ah, "k_OCP_HAP": k_oh}


def _ode(t, y, k_ao, k_ah, k_oh):
    """
    Three-state linear ODE.  y = [x_ACP, x_OCP, x_HAp], sum=1.
    Rates are first-order in the source phase.
    """
    x_ACP, x_OCP, x_HAp = y
    dx_ACP = -(k_ao + k_ah) * x_ACP
    dx_OCP =   k_ao * x_ACP  - k_oh * x_OCP
    dx_HAp =   k_ah * x_ACP  + k_oh * x_OCP
    return [dx_ACP, dx_OCP, dx_HAp]


def phase_evolution(
    t_days: np.ndarray,
    T_celsius: float = -20.0,
    pH: float = 7.81,
) -> dict:
    """
    Integrate ACP→OCP→HAp ODE from t=0 (pure ACP) over t_days.

    Returns dict with arrays x_ACP, x_OCP, x_HAp (fractions, sum=1).
    """
    kk = rate_constants(T_celsius, pH)
    sol = solve_ivp(
        _ode,
        t_span=(t_days[0], t_days[-1]),
        y0=[1.0, 0.0, 0.0],
        args=(kk["k_ACP_OCP"], kk["k_ACP_HAP"], kk["k_OCP_HAP"]),
        t_eval=t_days,
        method="Radau",
        rtol=1e-8,
        atol=1e-10,
    )
    return {
        "t_days": sol.t,
        "x_ACP":  sol.y[0],
        "x_OCP":  sol.y[1],
        "x_HAp":  sol.y[2],
        "rates":  kk,
    }


# ─────────────────────────────────────────────────────────────────────────────
# PART B: Post-thaw redissolution (Noyes-Whitney + Ostwald-Freundlich)
# ─────────────────────────────────────────────────────────────────────────────
#
# Ca recovery = fraction of precipitated Ca that re-dissolves by time t_thaw.
#
# For each mineral phase p:
#   dC/dt = D_p / (r_p + h) × (C_s,p(r_p) − C)
#
# where:
#   C_s,p(r) = C_s,bulk × exp(2γ_p Vm_p / (r R T))  [Ostwald-Freundlich]
#   h          = boundary layer thickness
#   r_p        = mean particle radius (function of storage time via LSW)
#   D_p        = diffusion coefficient of Ca²⁺ (corrected for temperature
#                and viscosity at 22°C after thaw)
#
# We compute Ca recovery = (C(t_thaw) − C_initial) / C_precipitated

# Physical constants
GAMMA = {          # surface tension (J/m²); Barralet 2004, Dorozhkin review
    "ACP":  0.07,
    "OCP":  0.15,
    "HAp":  0.20,
}
VM = {             # molar volume (m³/mol)
    "ACP":  67e-6,  # Ca₃(PO₄)₂, MW≈310, ρ≈2.5 g/cm³
    "OCP":  87e-6,  # Ca₄H(PO₄)₃, MW≈730, ρ≈2.6 g/cm³
    "HAp": 159e-6,  # Ca₅(PO₄)₃OH, MW≈502, ρ≈3.16 g/cm³
}
# Bulk solubility of each phase in mol-Ca/L at 22°C
C_SAT_BULK = {     # approximate; ACP ~ 0.4 mM, OCP ~ 0.05 mM, HAp ~ 1e-4 mM
    "ACP":  4.0e-4,    # mol/L = 0.4 mM  (Boskey 1973)
    "OCP":  5.0e-5,    # 0.05 mM  (Moreno 1960)
    "HAp":  1.0e-7,    # 0.0001 mM (McDowell 1977)
}

# Diffusion coefficient Ca²⁺ at 22°C in dilute water: 7.9e-10 m²/s
# In post-thaw solution (diluted glycerol, ~1% residual at 22°C):
#   D_eff ≈ 7.0e-10 m²/s (slight correction; Stokes-Einstein viscosity)
D_CA_22C = 7.0e-10   # m²/s

# Boundary layer thickness (m) per protocol
H_QUIESCENT  = 10e-6    # 10 µm quiescent; diffusion-limited for μm-scale particles
H_VORTEX     = 2e-6     # 2 µm vortex-mixed (literature: 1–5 µm)
H_COLD_REF   = 10e-6    # same geometry as quiescent but 4°C → slower D

# Diffusion coefficient correction for cold (2–8°C, use 4°C)
# Stokes-Einstein: D ∝ T/η; η_water(4°C)/η_water(22°C) ≈ 1.52/0.96 ≈ 1.58
D_CA_4C = D_CA_22C * (277.15 / 295.15) / 1.58   # ≈ 4.2e-10 m²/s

# LSW Ostwald ripening (r³ ∝ t):
#   <r>(t) = (r0³ + K_LSW × t)^(1/3)
# K_LSW estimated from Voorhees 1985 theory for Ca phosphate nanoparticles
# in concentrated solution at −20°C. Literature for ACP in solution at 25°C
# gives particle growth from ~30 nm to ~100 nm over ~24 h → K_LSW ≈ (100³-30³)/1 (nm³/h)
# At −20°C, viscosity correction reduces K_LSW by VISC_SCALE:
_K_LSW_ACP_25C = (100.0**3 - 30.0**3) / 24.0   # nm³/h at 25°C
K_LSW = {
    "ACP": _K_LSW_ACP_25C * VISC_SCALE,   # nm³/h at −20°C
    "OCP": _K_LSW_ACP_25C * VISC_SCALE * 0.3,   # OCP grows slower (higher γ)
    "HAp": _K_LSW_ACP_25C * VISC_SCALE * 0.1,   # HAp grows slowest
}
R0_NM = {"ACP": 30.0, "OCP": 50.0, "HAp": 80.0}   # initial mean radius (nm)


def mean_radius_nm(phase: str, t_storage_h: float) -> float:
    """Mean particle radius in nm after t_storage_h at −20°C via LSW theory."""
    r0 = R0_NM[phase]
    return (r0**3 + K_LSW[phase] * t_storage_h) ** (1.0 / 3.0)


def c_sat_ostwald(phase: str, r_nm: float, T_celsius: float = 22.0) -> float:
    """Effective solubility (mol/L) with Ostwald-Freundlich curvature correction."""
    r_m  = r_nm * 1e-9
    T_K  = T_celsius + 273.15
    expt = 2.0 * GAMMA[phase] * VM[phase] / (r_m * R_GAS * T_K)
    return C_SAT_BULK[phase] * np.exp(expt)


def ca_recovery_curve(
    storage_months: float,
    thaw_times_min: np.ndarray,
    protocol: str = "quiescent_22C",
    T_storage_C: float = -20.0,
    pH_storage: float = 7.81,
) -> np.ndarray:
    """
    Return Ca recovery fraction (0–1) at each time in thaw_times_min.

    storage_months — duration of storage at T_storage_C
    protocol       — "quiescent_22C", "vortex_22C", "cold_4C_48h"

    Ca recovery = fraction of precipitated Ca re-dissolved at thaw time t.

    The total Ca precipitated is computed from the phase fractions at
    the end of storage: each phase has a characteristic dissolution
    rate that drives recovery.
    """
    # Storage duration in days and hours
    t_storage_days = storage_months * 30.4375
    t_storage_h    = t_storage_days * 24.0

    # Phase fractions at end of storage
    sol = phase_evolution(
        np.array([0.0, t_storage_days]),
        T_celsius=T_storage_C,
        pH=pH_storage,
    )
    f_ACP = float(sol["x_ACP"][-1])
    f_OCP = float(sol["x_OCP"][-1])
    f_HAp = float(sol["x_HAp"][-1])

    # Protocol parameters
    if protocol == "quiescent_22C":
        D_eff = D_CA_22C
        h_m   = H_QUIESCENT
        T_thaw = 22.0
    elif protocol == "vortex_22C":
        D_eff = D_CA_22C
        h_m   = H_VORTEX
        T_thaw = 22.0
    elif protocol == "cold_4C_48h":
        D_eff = D_CA_4C
        h_m   = H_QUIESCENT
        T_thaw = 4.0
    else:
        raise ValueError(f"Unknown protocol: {protocol}")

    # For each phase: compute how much Ca re-dissolves by each thaw time.
    # Fraction dissolved from phase p: F_p(t) = 1 - exp(-lambda_p * t)
    # where lambda_p = D_eff / (r_p + h) / (V_particle) × A_s
    # Simplified: for nanoparticle (r < h), lambda ≈ D_eff / h / r × C_s / ΔC
    # Full Noyes-Whitney:  dF_p/dt = (D/h) * (C_s - C_solution) / C_precipitate
    # We integrate numerically for each phase using its dissolution rate.

    t_sec = thaw_times_min * 60.0   # convert to seconds
    recovery = np.zeros_like(thaw_times_min, dtype=float)

    for phase, f_phase in [("ACP", f_ACP), ("OCP", f_OCP), ("HAp", f_HAp)]:
        if f_phase < 1e-6:
            continue
        r_nm = mean_radius_nm(phase, t_storage_h)
        r_m  = r_nm * 1e-9
        c_s  = c_sat_ostwald(phase, r_nm, T_celsius=T_thaw)

        # Effective 1st-order dissolution rate constant (1/s):
        # lambda = D_eff * A/V / h_eff × (C_s / C_precipitate)
        # For a sphere: A/V = 3/r; h_eff = max(r, h_m) (diffusion layer)
        h_eff   = max(r_m, h_m)
        # C_precipitate ≈ 1/(VM) mol/m³ → mol/L = 1/(VM*1000)
        c_prec  = 1.0 / (VM[phase] * 1000.0)   # mol/L, molar concentration of solid
        lambda_p = D_eff * (3.0 / r_m) / h_eff * (c_s / c_prec)

        # Clamp to physically reasonable max dissolution rate
        lambda_p = min(lambda_p, 1.0)   # cap at 1/s (dissolves in seconds)

        # Dissolved fraction of this phase at time t: F_p = 1 - exp(-lambda_p * t)
        F_p = 1.0 - np.exp(-lambda_p * t_sec)
        recovery += f_phase * F_p

    # Clamp to [0, 1]
    return np.clip(recovery, 0.0, 1.0)


def ca_deficit_at_60min(storage_months: float, protocol: str = "quiescent_22C",
                        pH_storage: float = 7.81) -> float:
    """Ca deficit (fraction) at 60 min thaw for a given storage duration."""
    rec = ca_recovery_curve(
        storage_months,
        np.array([60.0]),
        protocol=protocol,
        pH_storage=pH_storage,
    )
    # Seeker measures Ca deficit = fraction of Ca NOT recovered by 60 min.
    # But total Ca in vial = serum Ca; precipitated fraction is a subset.
    # The deficit fraction as % of total Ca = f_precipitated × (1 - recovery)
    # f_precipitated ≈ (1 - alpha_Ca_free) × (total Ca).
    # At cryo state alpha_Ca = 0.117 → ~88% of Ca goes through the system.
    # But only the Ca that is actually trapped as slow-dissolving OCP/HAp
    # contributes to the 60-min deficit.
    f_precip = 0.88   # fraction of total Ca that precipitates at peak cryo
    return f_precip * (1.0 - float(rec[0]))


# ─────────────────────────────────────────────────────────────────────────────
# CSV output
# ─────────────────────────────────────────────────────────────────────────────

STORAGE_MONTHS = [1, 3, 6, 12, 24]
PROTOCOLS      = ["quiescent_22C", "vortex_22C", "cold_4C_48h"]
THAW_TIMES     = [10.0, 30.0, 60.0, 120.0, 240.0, 1440.0, 2880.0]  # minutes


def build_recovery_table(pH_storage: float = 7.81) -> list[dict]:
    rows = []
    for sm in STORAGE_MONTHS:
        sol = phase_evolution(
            np.array([0.0, sm * 30.4375]),
            pH=pH_storage,
        )
        dominant = "ACP" if sol["x_ACP"][-1] > 0.5 else \
                   "OCP" if sol["x_OCP"][-1] > sol["x_HAp"][-1] else "HAp"
        for prot in PROTOCOLS:
            tt = np.array(THAW_TIMES)
            rec = ca_recovery_curve(sm, tt, protocol=prot, pH_storage=pH_storage)
            row = {
                "storage_months":     sm,
                "thaw_protocol":      prot,
                "dominant_phase":     dominant,
                "x_ACP":              round(float(sol["x_ACP"][-1]), 4),
                "x_OCP":              round(float(sol["x_OCP"][-1]), 4),
                "x_HAp":              round(float(sol["x_HAp"][-1]), 4),
                "Ca_recovery_10min":  round(float(rec[0]), 4),
                "Ca_recovery_30min":  round(float(rec[1]), 4),
                "Ca_recovery_60min":  round(float(rec[2]), 4),
                "Ca_recovery_2h":     round(float(rec[3]), 4),
                "Ca_recovery_4h":     round(float(rec[4]), 4),
                "Ca_recovery_24h":    round(float(rec[5]), 4),
                "Ca_recovery_48h":    round(float(rec[6]), 4),
                "Ca_deficit_60min_%": round((1.0 - float(rec[2])) * 88.0, 2),
            }
            rows.append(row)
    return rows


def save_recovery_table(rows: list[dict], path: Path):
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point (kinetics + CSV; figures in separate functions below)
# ─────────────────────────────────────────────────────────────────────────────

def main_kinetics():
    print("Module 5: Ostwald ripening kinetics")
    print("=" * 60)

    kk = rate_constants(-20.0, 7.81)
    print(f"Rate constants at −20°C, pH 7.81:")
    print(f"  k(ACP→OCP) = {kk['k_ACP_OCP']:.4e} /day  "
          f"(t½ = {np.log(2)/kk['k_ACP_OCP']:.1f} d)")
    print(f"  k(ACP→HAp) = {kk['k_ACP_HAP']:.4e} /day")
    print(f"  k(OCP→HAp) = {kk['k_OCP_HAP']:.4e} /day  "
          f"(t½ = {np.log(2)/kk['k_OCP_HAP']:.1f} d)")
    print()

    # Phase fractions at key timepoints
    t_days = np.array([0, 7, 30, 90, 180, 365, 548, 730])
    sol = phase_evolution(t_days)
    print("Phase fractions vs storage time (pH 7.81, −20°C):")
    print(f"  {'Days':>6}  {'Mo':>4}  {'ACP':>6}  {'OCP':>6}  {'HAp':>6}  dominant")
    for i, td in enumerate(t_days):
        mo = td / 30.44
        dom = "ACP" if sol["x_ACP"][i] > 0.5 else \
              "OCP" if sol["x_OCP"][i] > sol["x_HAp"][i] else "HAp"
        print(f"  {td:>6.0f}  {mo:>4.1f}  "
              f"{sol['x_ACP'][i]:>6.3f}  "
              f"{sol['x_OCP'][i]:>6.3f}  "
              f"{sol['x_HAp'][i]:>6.3f}  {dom}")
    print()

    # Ca deficit at 60 min for key storage durations
    print("Ca deficit at 60-min thaw (quiescent 22°C):")
    for sm in [1, 3, 6, 12, 24]:
        deficit = ca_deficit_at_60min(sm)
        flag = " ← Seeker threshold" if abs(sm - 6) < 0.1 else ""
        print(f"  {sm:>2} months: {deficit*100:.1f}%{flag}")
    print()

    # Build and save CSV
    rows = build_recovery_table()
    save_recovery_table(rows, DATA_DIR / "module5_recovery_predictions.csv")
    return sol, rows


if __name__ == "__main__":
    main_kinetics()
