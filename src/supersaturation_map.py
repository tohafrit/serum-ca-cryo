"""
Module 4: Combined cryoconcentration × supersaturation map.

Sweeps cryoconcentration factor k (1–50) and pH (6.5–8.5) at T = −20°C,
computes SI for six mineral phases, and overlays the realistic cryo
trajectory under three CO₂-loss scenarios.  Produces fig03.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from pathlib import Path

from src.saturation_indices import (
    all_si,
    REFERENCE_COMPOSITION_mM,
)
from src.freezing_trajectory import (
    pH_trajectory,
    cryoconcentration_trajectory,
)

# ── Output directory ──────────────────────────────────────────────────────────

FIGURES_DIR = Path(__file__).parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# ── Grid parameters ───────────────────────────────────────────────────────────

# k on a log scale so low-k (near-physiological) region is well-resolved
K_GRID   = np.logspace(0, np.log10(50), 60)   # 1 → 50, 60 points
PH_GRID  = np.linspace(6.5, 8.5, 50)           # 50 pH steps

T_MAP    = -20.0          # all heatmaps at −20°C
I_BASE   = 0.1433         # mol/kg, baseline ionic strength (Module 2)

PHASES_PANEL = ["HAp", "OCP", "brushite", "ACP_loose", "ACP_tight", "calcite"]
PANEL_LABELS = {
    "HAp":      "HAp  Ca₅(PO₄)₃OH",
    "OCP":      "OCP  Ca₄H(PO₄)₃·2.5H₂O",
    "brushite": "Brushite  CaHPO₄·2H₂O",
    "ACP_loose":"ACP (loose)  Ca₃(PO₄)₂",
    "ACP_tight":"ACP (tight)  Ca₃(PO₄)₂",
    "calcite":  "Calcite  CaCO₃",
}

# ── Build composition dict for a given k ─────────────────────────────────────

BASE_COMP = dict(REFERENCE_COMPOSITION_mM)
# Protein in mM: 4.5 g/dL albumin (MW 66.5 kDa)
BASE_COMP["protein"] = 4.5 * 10.0 / 66.5 * 1000.0


def _comp_at_k(k: float) -> dict:
    """Return all concentrations scaled by k (simple cryoconcentration)."""
    return {sp: v * k for sp, v in BASE_COMP.items()}


# ── Compute SI grid ───────────────────────────────────────────────────────────

def compute_si_grid() -> dict:
    """
    Returns dict[phase] → 2-D array shape (len(K_GRID), len(PH_GRID)).
    Rows = k, columns = pH.
    """
    nk, npH = len(K_GRID), len(PH_GRID)
    grids = {ph: np.empty((nk, npH)) for ph in PHASES_PANEL}

    for i, k in enumerate(K_GRID):
        comp = _comp_at_k(k)
        I    = I_BASE * k
        for j, pH in enumerate(PH_GRID):
            si = all_si(comp, pH=pH, T_celsius=T_MAP, I=I)
            for ph in PHASES_PANEL:
                grids[ph][i, j] = si[ph]

    return grids


# ── Cryo trajectory overlay ───────────────────────────────────────────────────

# CO₂ loss fractions and their display colours
CO2_SCENARIOS = [
    (0.0,  "#1f77b4", "0% CO₂ loss",  "-"),
    (0.5,  "#ff7f0e", "50% CO₂ loss", "--"),
    (0.9,  "#d62728", "90% CO₂ loss", ":"),
]


def _trajectory_k_pH(co2_loss_frac: float):
    """
    Return (k_arr, pH_arr) for the cryo trajectory at T = −20°C as k varies.

    k is taken directly from K_GRID; pH is computed via pH_trajectory which
    uses the H-H equation with T- and I-corrected pKa₁.
    """
    df = pH_trajectory(
        K_GRID,
        BASE_COMP,
        co2_loss_fractions=(co2_loss_frac,),
        T_celsius=T_MAP,
    )
    col = f"pH_co2loss_{int(round(co2_loss_frac*100)):02d}pct"
    return K_GRID, df[col].values


# ── Figure ────────────────────────────────────────────────────────────────────

# Symmetric diverging limits per panel (SI=0 is saturation boundary)
VLIM = {
    "HAp":       (-2, 14),
    "OCP":       (-10, 8),
    "brushite":  (-3, 3),
    "ACP_loose": (-8, 4),
    "ACP_tight": (-5, 5),
    "calcite":   (-3, 3),
}

# Use RdBu_r: blue=undersaturated, red=supersaturated
CMAP = "RdBu_r"


def plot_si_heatmaps(grids: dict, out_path: Path):
    """
    2×3 panel heatmap figure.

    x-axis: pH (6.5 – 8.5)
    y-axis: k (log scale, 1 – 50)
    colour: SI (white = SI=0 saturation boundary)
    overlay: cryo trajectory for 0/50/90% CO₂ loss
    """
    fig, axes = plt.subplots(
        2, 3,
        figsize=(7.09, 5.0),   # 180 mm × 127 mm
        constrained_layout=True,
        sharex=True, sharey=True,
    )

    ph_ext  = [PH_GRID[0], PH_GRID[-1]]
    k_ext   = [K_GRID[0],  K_GRID[-1]]
    # extent for imshow: [left, right, bottom, top] in data coords
    # imshow treats row 0 as top; we want k increasing upward → flip
    extent = [ph_ext[0], ph_ext[1], np.log10(k_ext[0]), np.log10(k_ext[1])]

    # Pre-compute trajectory arrays once
    traj = {frac: _trajectory_k_pH(frac) for frac, *_ in CO2_SCENARIOS}

    for ax, phase in zip(axes.flat, PHASES_PANEL):
        Z    = grids[phase]          # shape (nk, npH), k=rows, pH=cols
        vmin, vmax = VLIM[phase]

        # Diverging normalisation centred at SI=0
        norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

        # imshow: rows → y-axis; we flip Z vertically so small k is at bottom
        im = ax.imshow(
            Z[::-1, :],
            aspect="auto",
            extent=extent,
            norm=norm,
            cmap=CMAP,
            origin="upper",   # after flip, first row of flipped Z = largest k
            interpolation="bilinear",
        )

        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        cb.set_label("SI", fontsize=7)
        cb.ax.tick_params(labelsize=6)

        # SI = 0 contour (saturation boundary)
        # Need to work in log10(k) coordinates for the y-axis
        log10_k = np.log10(K_GRID)
        ax.contour(
            PH_GRID, log10_k, Z,
            levels=[0.0],
            colors="black",
            linewidths=1.2,
        )

        # Cryo trajectory overlays
        for frac, colour, label, ls in CO2_SCENARIOS:
            k_arr, pH_arr = traj[frac]
            ax.plot(
                pH_arr,
                np.log10(k_arr),
                color=colour, ls=ls, lw=1.4,
                zorder=5,
            )

        # Mark physiological starting point
        ax.plot(7.4, np.log10(1.0), "w*", ms=7, zorder=6)

        # Axis cosmetics
        ax.set_title(PANEL_LABELS[phase], fontsize=7.5, pad=3)
        ax.set_xlim(PH_GRID[0], PH_GRID[-1])
        ax.set_ylim(np.log10(K_GRID[0]), np.log10(K_GRID[-1]))

        # y-ticks in natural k units
        k_ticks = [1, 2, 5, 10, 20, 50]
        ax.set_yticks([np.log10(k) for k in k_ticks])
        ax.set_yticklabels([str(k) for k in k_ticks], fontsize=6)
        ax.tick_params(axis="x", labelsize=6)

    # Shared axis labels
    for ax in axes[1, :]:
        ax.set_xlabel("pH", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("Conc. factor k", fontsize=8)

    # Shared trajectory legend — placed outside the grid at bottom
    legend_handles = [
        Line2D([0], [0], color=col, ls=ls, lw=1.4, label=lbl)
        for _, col, lbl, ls in CO2_SCENARIOS
    ] + [
        Line2D([0], [0], color="black", lw=1.2, label="SI = 0 (saturation)"),
        Line2D([0], [0], marker="*", color="white", markerfacecolor="white",
               markersize=7, ls="None", label="Physiological start"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        fontsize=6.5,
        frameon=True,
        bbox_to_anchor=(0.5, -0.06),
    )

    fig.suptitle(
        "Saturation indices at −20 °C: cryoconcentration factor k vs pH",
        fontsize=9, y=1.01,
    )

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    print("Module 4: Supersaturation map")
    print("=" * 60)
    print(f"Grid: {len(K_GRID)} k-values × {len(PH_GRID)} pH-values = "
          f"{len(K_GRID)*len(PH_GRID):,} SI evaluations")
    print(f"Temperature: {T_MAP} °C")
    print()

    print("Computing SI grid …")
    grids = compute_si_grid()

    # Quick summary: SI(HAp) at physiological and cryo points
    idx_k1   = 0
    idx_k558 = int(np.argmin(np.abs(K_GRID - 5.58)))
    idx_pH74 = int(np.argmin(np.abs(PH_GRID - 7.4)))
    idx_pH78 = int(np.argmin(np.abs(PH_GRID - 7.81)))

    si_hap_phys = grids["HAp"][idx_k1,   idx_pH74]
    si_hap_cryo = grids["HAp"][idx_k558, idx_pH78]

    print(f"  SI(HAp) at k=1.0,  pH=7.4  → {si_hap_phys:+.2f}")
    print(f"  SI(HAp) at k=5.58, pH=7.81 → {si_hap_cryo:+.2f}")
    print()

    out = FIGURES_DIR / "fig03_si_heatmaps.png"
    print("Generating fig03 …")
    plot_si_heatmaps(grids, out)
    print()
    print("Done.")
    return grids


if __name__ == "__main__":
    main()
