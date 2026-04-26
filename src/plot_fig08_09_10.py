"""fig08, fig09, fig10: Intervention efficacy figures."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

from src.interventions import (
    run_all_scenarios, SCENARIOS, scenario_stats, STORAGE_MONTHS, THRESHOLD,
)

FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Display names and colours
DISPLAY = {
    "baseline":           ("Baseline (pH 8.0)",           "#555555"),
    "+degas_10":          ("+Degas (10% CO₂)",            "#2196F3"),
    "+crf_2C":            ("+CRF 2°C/min",                "#4CAF50"),
    "+vortex_30s":        ("+Vortex 30 s",                "#FF9800"),
    "+combined":          ("+All three (combined)",       "#d62728"),
    "+combined_plus":     ("+Combined+ (double vortex)",  "#9C27B0"),
    "+seeker_workaround": ("Seeker 48-h cold soak*",      "#795548"),
}
# 5 standard-thaw scenarios for fig08 histograms (comparable 60-min thaw)
SHOW_SCENARIOS    = ["baseline", "+degas_10", "+crf_2C", "+vortex_30s", "+combined"]
# Extended set for fig09 bar chart (includes the two new scenarios)
SHOW_SCENARIOS_09 = list(DISPLAY.keys())


# ── fig08: deficit distribution histograms at 12 months ──────────────────────

def plot_fig08(all_results: dict):
    fig, axes = plt.subplots(1, 5, figsize=(13.0, 3.5),
                             sharey=False, constrained_layout=True)

    j12 = STORAGE_MONTHS.index(12)
    for ax, name in zip(axes, SHOW_SCENARIOS):
        label, col = DISPLAY[name]
        d = all_results[name][:, j12] * 100.0
        frac = (d > THRESHOLD * 100).mean()
        ax.hist(d, bins=50, color=col, alpha=0.80, edgecolor="none")
        ax.axvline(THRESHOLD * 100, color="crimson", lw=1.5, ls="--")
        ax.set_title(label, fontsize=7.5, pad=4)
        ax.set_xlabel("Ca deficit (%)", fontsize=7.5)
        ax.tick_params(labelsize=7)
        ax.text(0.97, 0.97, f"{frac:.0%}\n>4%",
                transform=ax.transAxes, fontsize=9,
                ha="right", va="top", color="crimson", fontweight="bold")

    axes[0].set_ylabel("Vials", fontsize=8)
    fig.suptitle("Post-thaw Ca deficit distribution at 12-month storage\n"
                 "N = 10,000 vials per scenario, 22°C quiescent 60-min thaw",
                 fontsize=9)
    out = FIGURES_DIR / "fig08_intervention_distributions.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── fig09: efficacy bar chart ─────────────────────────────────────────────────

def plot_fig09(rows: list[dict]):
    rows5 = [r for r in rows if r["scenario"] in SHOW_SCENARIOS_09]
    names  = [DISPLAY[r["scenario"]][0] for r in rows5]
    cols   = [DISPLAY[r["scenario"]][1] for r in rows5]
    f6     = [r["frac_above_4pct_6mo"]  * 100 for r in rows5]
    f12    = [r["frac_above_4pct_12mo"] * 100 for r in rows5]

    x   = np.arange(len(names))
    w   = 0.30
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12.0, 4.2),
                                             constrained_layout=True)

    # ── Left: fraction above threshold ──────────────────────────────────────
    bars6  = ax_left.bar(x - w/2, f6,  w, color=cols, alpha=0.9, label="6 months")
    bars12 = ax_left.bar(x + w/2, f12, w, color=cols, alpha=0.55, label="12 months",
                         hatch="///", edgecolor="white", linewidth=0.5)
    ax_left.set_ylabel("Vials with Ca deficit > 4%  (%)", fontsize=8.5)
    ax_left.set_xticks(x)
    ax_left.set_xticklabels(names, rotation=30, ha="right", fontsize=7)
    ax_left.axhline(100, color="gray", lw=0.5, ls=":")
    ax_left.set_ylim(0, 112)
    ax_left.set_title("Fraction of vials above 4% threshold", fontsize=9)
    # Value labels — show all, "0%" explicitly for zero bars
    for bar in list(bars6) + list(bars12):
        h = bar.get_height()
        if h > 0.5:
            ax_left.text(bar.get_x() + bar.get_width()/2, h + 1,
                         f"{h:.0f}%", ha="center", va="bottom", fontsize=6.0)
        else:
            ax_left.text(bar.get_x() + bar.get_width()/2, 1.5,
                         "0%", ha="center", va="bottom", fontsize=6.0,
                         color="gray")
    ax_left.legend(fontsize=7.5, loc="upper right")
    ax_left.tick_params(labelsize=7)
    ax_left.annotate("* 48-h cold soak;\ndifferent thaw protocol",
                     xy=(0.98, 0.02), xycoords="axes fraction",
                     fontsize=6, ha="right", va="bottom", color="#795548",
                     style="italic")

    # ── Right: mean deficit reduction vs baseline ────────────────────────────
    baseline_12 = rows5[0]["mean_deficit_12mo_pct"]
    reductions  = [baseline_12 - r["mean_deficit_12mo_pct"] for r in rows5]
    bar2 = ax_right.bar(x, reductions, color=cols, alpha=0.85)
    ax_right.set_ylabel("Mean deficit reduction at 12 months  (pp)", fontsize=8.5)
    ax_right.set_xticks(x)
    ax_right.set_xticklabels(names, rotation=30, ha="right", fontsize=7)
    ax_right.set_title("Mean deficit reduction vs baseline (12 months)", fontsize=9)
    ax_right.axhline(0, color="black", lw=0.5)
    for bar in bar2:
        h = bar.get_height()
        ax_right.text(bar.get_x() + bar.get_width()/2, h + 0.05,
                      f"{h:.1f}", ha="center", va="bottom", fontsize=7.0)
    ax_right.tick_params(labelsize=7)

    out = FIGURES_DIR / "fig09_intervention_efficacy.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


# ── fig10: Pareto plot (efficacy vs implementation complexity) ────────────────

COMPLEXITY = {   # 1=very easy, 5=very hard (capex + operational)
    "baseline":           (1, "Baseline"),
    "+degas_50":          (2, "+Degas 50%"),
    "+degas_10":          (2, "+Degas 10%"),
    "+crf_1C":            (4, "+CRF 1°C/min"),
    "+crf_2C":            (5, "+CRF 2°C/min"),
    "+vortex_30s":        (2, "+Vortex 30s"),
    "+vortex_60s":        (2, "+Vortex 60s"),
    "+combined":          (5, "+Combined"),
    "+combined_plus":     (5, "+Combined+"),
    "+seeker_workaround": (1, "Seeker 48-h*"),
}

def plot_fig10(rows: list[dict]):
    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    for row in rows:
        name = row["scenario"]
        comp, label = COMPLEXITY[name]
        f6   = row["frac_above_4pct_6mo"]  * 100
        col  = DISPLAY.get(name, (None, "#888888"))[1]
        size = 120 if name in DISPLAY else 60

        ax.scatter(comp, f6, s=size, color=col, zorder=5,
                   edgecolors="white", linewidths=0.7)
        offset_x = 0.12 if comp < 4 else -0.12
        ha        = "left" if comp < 4 else "right"
        ax.annotate(label, (comp, f6),
                    textcoords="offset points", xytext=(5, 3),
                    fontsize=7, ha=ha, color="#222222")

    ax.axhline(10, color="gray", ls=":", lw=1.0)
    ax.text(0.5, 10.5, "10% target", fontsize=7, color="gray")
    ax.set_xlabel("Implementation complexity (1=trivial, 5=capital equipment)", fontsize=9)
    ax.set_ylabel("Fraction of vials >4% at 6 months  (%)", fontsize=9)
    ax.set_title("Pareto: efficacy vs implementation cost", fontsize=9)
    ax.set_xlim(0.5, 5.8)
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xticklabels(["1\nSOP change", "2\nMinor capex",
                         "3", "4\nMajor capex", "5\nFull system"], fontsize=7)
    ax.tick_params(labelsize=7)
    ax.set_ylim(-3, 65)

    # Highlight Pareto frontier (manually — scenarios not dominated)
    # Vortex 60s: complexity 2, efficacy 6.3% — not dominated
    # Combined: complexity 5, efficacy 0.4% — not dominated
    ax.annotate("← Best single\n   intervention",
                xy=(2, 6.3), xytext=(2.5, 15),
                fontsize=7, color="#FF9800",
                arrowprops=dict(arrowstyle="->", color="#FF9800", lw=1.0))

    fig.tight_layout()
    out = FIGURES_DIR / "fig10_pareto.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Running all intervention scenarios …")
    all_results = run_all_scenarios(seed=42)
    rows = scenario_stats(all_results)

    print("\nGenerating figures …")
    plot_fig08(all_results)
    plot_fig09(rows)
    plot_fig10(rows)
    print("Done.")
