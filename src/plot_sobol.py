"""
Sobol sensitivity analysis for Module 6.

Uses SALib Saltelli sampling. Each parameter is sampled uniformly over
its plausible range (log-uniform for lognormal variables). Computes
first-order (S1) and total-order (ST) Sobol indices for Ca deficit
at 60-min thaw at 12-month storage.
"""

import numpy as np
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from SALib.sample import saltelli
from SALib.analyze import sobol as sobol_analyze

from src.vial_simulation import compute_vial_deficit, STORAGE_MONTHS

FIGURES_DIR = Path(__file__).parent.parent / "figures"
DATA_DIR    = Path(__file__).parent.parent / "data"

STORAGE_MONTH = 12    # evaluate at 12-month storage

# ── Sobol problem definition ──────────────────────────────────────────────────
# Parameters varied uniformly within physically plausible bounds.
# Lognormal parameters (nucleation_delay, particle_size, etc.) are sampled
# in log₁₀ space and transformed back before calling compute_vial_deficit.

PROBLEM = {
    "num_vars": 9,
    "names": [
        "local_k_factor",          # spatial freeze heterogeneity
        "log10_nuc_delay_days",    # log10(nucleation_delay_days)
        "log10_particle_nm",       # log10(initial ACP particle size, nm)
        "storage_T_C",             # actual freezer temperature
        "albumin_g_dL",            # serum albumin (donor variability)
        "cryo_purity",             # glycerol lot purity
        "log10_glass_density",     # log10(glass surface site density)
        "fill_volume_mL",          # vial fill volume
        "log10_freezing_rate",     # log10(freezing rate, °C/min)
    ],
    "bounds": [
        [0.6,  1.6],              # local_k_factor
        [np.log10(1.0), np.log10(900.0)],   # log10(nuc_delay_days): 1–900 d
        [np.log10(10.0), np.log10(300.0)],  # log10(particle_nm): 10–300 nm
        [-24.0, -16.0],           # storage_T_C
        [3.5,   7.5],             # albumin_g_dL
        [0.94,  1.06],            # cryo_purity
        [np.log10(0.05), np.log10(5.0)],    # log10(glass_density)
        [6.5,   8.5],             # fill_volume_mL
        [np.log10(0.3), np.log10(5.0)],     # log10(freezing_rate)
    ],
}

PARAM_LABELS = {
    "local_k_factor":       "Local cryo-factor k",
    "log10_nuc_delay_days": "Nucleation delay (days)",
    "log10_particle_nm":    "Initial ACP particle size",
    "storage_T_C":          "Storage temperature (°C)",
    "albumin_g_dL":         "Serum albumin (g/dL)",
    "cryo_purity":          "Cryoprotectant purity",
    "log10_glass_density":  "Glass surface site density",
    "fill_volume_mL":       "Fill volume (mL)",
    "log10_freezing_rate":  "Freezing rate (°C/min)",
}


def _evaluate(X: np.ndarray) -> np.ndarray:
    """Evaluate Ca deficit for each row in the Saltelli sample matrix X."""
    Y = np.empty(len(X))
    for i, row in enumerate(X):
        (local_k, log_nuc, log_part, T_C,
         albumin, cryo, log_glass, fill, log_freeze) = row
        Y[i] = compute_vial_deficit(
            storage_months        = STORAGE_MONTH,
            local_k_factor        = local_k,
            nucleation_delay_days = 10.0 ** log_nuc,
            particle_size_nm      = 10.0 ** log_part,
            storage_T_C           = T_C,
            albumin_g_dL          = albumin,
            cryo_purity           = cryo,
            glass_density         = 10.0 ** log_glass,
            fill_volume_mL        = fill,
            freezing_rate         = 10.0 ** log_freeze,
        )
    return Y


def run_sobol(n_base: int = 1024) -> dict:
    """
    Run Sobol analysis with n_base × (2D+2) = n_base × 20 model evaluations.
    Returns SALib analysis result dict.
    """
    print(f"  Sobol sampling: N_base={n_base}, "
          f"evaluations={n_base * (2*PROBLEM['num_vars']+2):,}")
    X = saltelli.sample(PROBLEM, n_base, calc_second_order=False)
    print(f"  Evaluating {len(X):,} samples …")
    Y = _evaluate(X)
    Si = sobol_analyze.analyze(PROBLEM, Y, calc_second_order=False, print_to_console=False)
    return Si


def save_sobol_csv(Si: dict, path: Path):
    names = PROBLEM["names"]
    rows  = []
    for i, name in enumerate(names):
        rows.append({
            "parameter": name,
            "label":     PARAM_LABELS[name],
            "S1":        round(float(Si["S1"][i]), 4),
            "S1_conf":   round(float(Si["S1_conf"][i]), 4),
            "ST":        round(float(Si["ST"][i]), 4),
            "ST_conf":   round(float(Si["ST_conf"][i]), 4),
        })
    # Sort by ST descending
    rows.sort(key=lambda r: r["ST"], reverse=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved: {path}")
    return rows


def plot_sobol(rows: list[dict]):
    labels = [PARAM_LABELS[r["parameter"]] for r in rows]
    S1     = [r["S1"]    for r in rows]
    ST     = [r["ST"]    for r in rows]
    S1_ci  = [r["S1_conf"] for r in rows]
    ST_ci  = [r["ST_conf"] for r in rows]

    y   = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(6.5, 4.0))

    ax.barh(y + 0.18, ST, height=0.35, color="#d62728", alpha=0.8,
            label="Total-order $S_T$", xerr=ST_ci, capsize=3)
    ax.barh(y - 0.18, S1, height=0.35, color="#1f77b4", alpha=0.8,
            label="First-order $S_1$",  xerr=S1_ci, capsize=3)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Sobol index", fontsize=9)
    ax.set_title(f"Sobol sensitivity indices — Ca deficit at {STORAGE_MONTH}-month storage",
                 fontsize=9)
    ax.legend(fontsize=8, loc="lower right")
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlim(-0.05, 1.05)
    ax.tick_params(labelsize=7)

    fig.tight_layout()
    out = FIGURES_DIR / "figS4_sobol.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Sobol sensitivity analysis")
    print("=" * 50)
    Si   = run_sobol(n_base=1024)
    rows = save_sobol_csv(Si, DATA_DIR / "module6_sobol_indices.csv")

    print("\nTop parameters by ST:")
    for r in rows[:5]:
        print(f"  {r['label']:<35}  S1={r['S1']:.3f}  ST={r['ST']:.3f}")

    plot_sobol(rows)
    print("\nDone.")
