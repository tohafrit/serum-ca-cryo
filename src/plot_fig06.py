"""fig06: Vial deficit distributions at 4 storage durations."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from src.vial_simulation import run_monte_carlo, STORAGE_MONTHS, THRESHOLD

FIGURES_DIR = Path(__file__).parent.parent / "figures"

SHOW_MONTHS = [3, 6, 12, 24]
COLORS      = ["#2ca02c", "#ff7f0e", "#d62728", "#9467bd"]


def plot_fig06(result: dict):
    deficits = result["deficits"]
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.0), constrained_layout=True)

    for ax, sm, col in zip(axes.flat, SHOW_MONTHS, COLORS):
        j  = STORAGE_MONTHS.index(sm)
        d  = deficits[:, j] * 100.0   # → percent

        frac_above = (d > THRESHOLD * 100).mean()
        mean_d     = d.mean()
        med_d      = float(np.median(d))

        ax.hist(d, bins=60, color=col, alpha=0.75, edgecolor="none")
        ax.axvline(THRESHOLD * 100, color="crimson", lw=1.8, ls="--")
        ax.text(THRESHOLD * 100 + 0.1, ax.get_ylim()[1] * 0.5,
                f"detection ~0.5%\n{frac_above:.0%} affected",
                fontsize=7.5, color="crimson", va="center")

        ax.set_title(f"{sm}-month storage", fontsize=9)
        ax.set_xlabel("Ca deficit at 60-min thaw (%)", fontsize=8)
        ax.set_ylabel("Vials", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.text(0.97, 0.95,
                f"mean {mean_d:.1f}%\nmedian {med_d:.1f}%",
                transform=ax.transAxes, fontsize=7.5,
                ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    fig.suptitle(
        "Vial-to-vial distribution of post-thaw Ca deficit\n"
        f"N = 10,000 vials (50 batches × 200 vials), 22°C quiescent 60-min thaw",
        fontsize=9,
    )
    out = FIGURES_DIR / "fig06_vial_distributions.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    result, _ = __import__("src.vial_simulation", fromlist=["main_mc"]).main_mc()
    plot_fig06(result)
