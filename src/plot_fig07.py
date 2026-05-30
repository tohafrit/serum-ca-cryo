"""fig07: Batch-to-batch variability violin plot at 12 months."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from src.vial_simulation import (
    run_monte_carlo, STORAGE_MONTHS, THRESHOLD,
    N_BATCHES, N_PER_BATCH,
)

FIGURES_DIR = Path(__file__).parent.parent / "figures"


def plot_fig07(result: dict):
    deficits     = result["deficits"]
    batch_params = result["batch_params"]   # [albumin, cryo_purity, freezing_rate, storage_T, nuc_batch_scale]
    vial_params  = result["vial_params"]    # [batch_idx, ...]

    j12 = STORAGE_MONTHS.index(12)

    # Collect per-batch deficit arrays
    batch_deficits = []
    for b in range(N_BATCHES):
        mask = vial_params[:, 0].astype(int) == b
        batch_deficits.append(deficits[mask, j12] * 100.0)

    # Sort batches by median deficit (ascending) for visual clarity
    medians    = [np.median(d) for d in batch_deficits]
    sort_order = np.argsort(medians)
    sorted_d   = [batch_deficits[i] for i in sort_order]
    sorted_nuc = batch_params[sort_order, 4]   # nucleation_batch_scale (lower → shorter delay → worse)

    fig, (ax_main, ax_ann) = plt.subplots(
        1, 2, figsize=(9.0, 4.5),
        gridspec_kw={"width_ratios": [3, 1]},
        constrained_layout=True,
    )

    # ── Main violin ──────────────────────────────────────────────────────────
    parts = ax_main.violinplot(
        sorted_d,
        positions=np.arange(N_BATCHES),
        widths=0.7,
        showmedians=True,
        showextrema=False,
    )
    for pc in parts["bodies"]:
        pc.set_facecolor("#5B8DB8")
        pc.set_alpha(0.65)
    parts["cmedians"].set_color("#1f2b3a")
    parts["cmedians"].set_linewidth(1.2)

    ax_main.axhline(THRESHOLD * 100, color="crimson", lw=1.4, ls="--")
    ax_main.text(-0.5, THRESHOLD * 100 + 0.1, "4% threshold",
                 fontsize=7.5, color="crimson")

    ax_main.set_xlabel("Batch (sorted by median deficit)", fontsize=9)
    ax_main.set_ylabel("Ca deficit at 60-min thaw (%)", fontsize=9)
    ax_main.set_title("Batch-to-batch variability — 12-month storage\n"
                       "50 batches × 200 vials each", fontsize=9)
    ax_main.set_xlim(-1, N_BATCHES)
    ax_main.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_main.tick_params(labelsize=7)

    # ── Annotation panel: batch nucleation scale vs median deficit ───────────
    ax_ann.scatter(
        sorted_nuc,
        [np.median(d) for d in sorted_d],
        c=[np.median(d) for d in sorted_d],
        cmap="RdYlGn_r",
        s=30, alpha=0.8, edgecolors="none",
    )
    ax_ann.set_xlabel("Batch nuc. scale\n(low = faster nucleation)", fontsize=8)
    ax_ann.set_ylabel("Batch median deficit (%)", fontsize=8)
    ax_ann.set_title("Nucleation scale\nvs deficit", fontsize=8)
    ax_ann.axhline(THRESHOLD * 100, color="crimson", lw=1.0, ls="--")
    ax_ann.tick_params(labelsize=7)

    out = FIGURES_DIR / "fig07_batch_variability.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  Saved: {out}")

    # Print batch statistics
    all_medians = [np.median(d) for d in batch_deficits]
    print(f"  Batch median deficit at 12 months: "
          f"min {min(all_medians):.1f}%, "
          f"max {max(all_medians):.1f}%, "
          f"IQR [{np.percentile(all_medians,25):.1f}–{np.percentile(all_medians,75):.1f}]%")


if __name__ == "__main__":
    from src.vial_simulation import main_mc
    result, _ = main_mc()
    plot_fig07(result)
