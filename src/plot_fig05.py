"""fig05: Post-thaw Ca recovery curves — 2×2 panel."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from src.ripening_kinetics import ca_recovery_curve

FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Axis: time in minutes up to 60 min; then extended for 48 h panel
T_FAST = np.linspace(0, 60, 300)      # 0–60 min
T_SLOW = np.linspace(0, 2880, 1000)   # 0–48 h

COLORS = {1: "#2ca02c", 3: "#1f77b4", 6: "#ff7f0e", 12: "#d62728"}
MONTHS = [1, 3, 6, 12]

F_PRECIP = 0.88   # fraction of total Ca that precipitates


def _deficit_label(storage_mo, t_min, protocol="quiescent_22C"):
    """Ca deficit % at t_min for a given storage duration."""
    rec = ca_recovery_curve(storage_mo, np.array([t_min]), protocol=protocol)
    return F_PRECIP * (1.0 - float(rec[0])) * 100.0


def plot_fig05():
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.5), constrained_layout=True)

    # ── Top-left: quiescent 22°C, 0–60 min, all storage durations ──────────
    ax = axes[0, 0]
    for sm in MONTHS:
        rec = ca_recovery_curve(sm, T_FAST, protocol="quiescent_22C")
        deficit_pct = F_PRECIP * (1.0 - rec) * 100.0
        ax.plot(T_FAST, deficit_pct, color=COLORS[sm], lw=1.8,
                label=f"{sm} month{'s' if sm>1 else ''}")
    ax.axhline(4.0, color="crimson", ls="--", lw=2.0)
    ax.text(1, 4.3, "4% Seeker threshold", fontsize=7, color="crimson", fontweight="bold")
    ax.set_title("22°C quiescent (standard protocol)", fontsize=8.5)
    ax.set_xlabel("Time after thaw (min)", fontsize=8)
    ax.set_ylabel("Ca deficit (% of total)", fontsize=8)
    ax.set_xlim(0, 60); ax.set_ylim(0, None)
    ax.legend(fontsize=7.5, title="Storage", title_fontsize=7.5)
    ax.tick_params(labelsize=7)

    # ── Top-right: 22°C with 30-s vortex at t=30 min ───────────────────────
    ax = axes[0, 1]
    t_vortex = T_FAST.copy()
    for sm in MONTHS:
        rec_q = ca_recovery_curve(sm, t_vortex, protocol="quiescent_22C")
        rec_v = ca_recovery_curve(sm, t_vortex, protocol="vortex_22C")
        # Vortex applied at t=30 min: before 30 min use quiescent, after use vortex
        rec_combined = np.where(t_vortex <= 30.0, rec_q, rec_v)
        deficit_pct  = F_PRECIP * (1.0 - rec_combined) * 100.0
        ax.plot(t_vortex, deficit_pct, color=COLORS[sm], lw=1.8,
                label=f"{sm} month{'s' if sm>1 else ''}")
    ax.axvline(30.0, color="black", ls="--", lw=1.0)
    ax.text(31, ax.get_ylim()[1] * 0.05 if ax.get_ylim()[1] > 0 else 1,
            "vortex", fontsize=7)
    ax.axhline(4.0, color="crimson", ls="--", lw=2.0)
    ax.set_title("22°C + 30-s vortex at 30 min", fontsize=8.5)
    ax.set_xlabel("Time after thaw (min)", fontsize=8)
    ax.set_ylabel("Ca deficit (% of total)", fontsize=8)
    ax.set_xlim(0, 60); ax.set_ylim(0, None)
    ax.legend(fontsize=7.5, title="Storage", title_fontsize=7.5)
    ax.tick_params(labelsize=7)

    # ── Bottom-left: 2–8°C quiescent, 0–48 h ───────────────────────────────
    ax = axes[1, 0]
    T_h = T_SLOW / 60.0   # in hours for x-axis
    for sm in MONTHS:
        rec = ca_recovery_curve(sm, T_SLOW, protocol="cold_4C_48h")
        deficit_pct = F_PRECIP * (1.0 - rec) * 100.0
        ax.plot(T_h, deficit_pct, color=COLORS[sm], lw=1.8,
                label=f"{sm} month{'s' if sm>1 else ''}")
    ax.axhline(4.0, color="crimson", ls="--", lw=2.0)
    ax.annotate("60-min snapshot in cold protocol:\nkinetics still in progress",
                xy=(1.0, _deficit_label(12, 60.0, "cold_4C_48h")),
                xytext=(12, 10), fontsize=6.5, color="#555555",
                arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8))
    ax.set_title("2–8°C quiescent (48-h equilibration)", fontsize=8.5)
    ax.set_xlabel("Time after thaw (h)", fontsize=8)
    ax.set_ylabel("Ca deficit (% of total)", fontsize=8)
    ax.set_xlim(0, 48); ax.set_ylim(0, None)
    ax.legend(fontsize=7.5, title="Storage", title_fontsize=7.5)
    ax.tick_params(labelsize=7)

    # ── Bottom-right: 12-month sample, all protocols compared at 0–60 min ──
    ax = axes[1, 1]
    t_ext = np.linspace(0, 60, 300)
    for prot, lbl, col in [
        ("quiescent_22C", "22°C quiescent", "#d62728"),
        ("vortex_22C",    "22°C + vortex",  "#ff7f0e"),
        ("cold_4C_48h",   "2–8°C 48 h",     "#1f77b4"),
    ]:
        rec = ca_recovery_curve(12, t_ext, protocol=prot)
        deficit_pct = F_PRECIP * (1.0 - rec) * 100.0
        ax.plot(t_ext, deficit_pct, color=col, lw=1.8, label=lbl)
    ax.axhline(4.0, color="crimson", ls="--", lw=2.0)
    ax.text(1, 4.3, "4%", fontsize=7, color="crimson", fontweight="bold")
    ax.set_title("12-month sample: protocol comparison", fontsize=8.5)
    ax.set_xlabel("Time after thaw (min)", fontsize=8)
    ax.set_ylabel("Ca deficit (% of total)", fontsize=8)
    ax.set_xlim(0, 60); ax.set_ylim(0, None)
    ax.legend(fontsize=7.5)
    ax.tick_params(labelsize=7)

    fig.suptitle("Post-thaw Ca recovery — Noyes-Whitney / Ostwald-Freundlich model",
                 fontsize=9.5)
    out = FIGURES_DIR / "fig05_redissolution.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  Saved: {out}")

    # Print key numbers for verification
    print("\nKey deficits at 60 min (quiescent 22°C):")
    for sm in MONTHS:
        d = _deficit_label(sm, 60.0)
        print(f"  {sm:>2} months: {d:.1f}%")
    print("\nDeficits at 60 min, 12-month sample:")
    for prot in ["quiescent_22C", "vortex_22C", "cold_4C_48h"]:
        d = _deficit_label(12, 60.0, protocol=prot)
        print(f"  {prot}: {d:.1f}%")


if __name__ == "__main__":
    plot_fig05()
