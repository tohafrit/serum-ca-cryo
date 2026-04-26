"""fig04: Phase evolution during storage at −20°C (stacked area chart)."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from src.ripening_kinetics import phase_evolution

FIGURES_DIR = Path(__file__).parent.parent / "figures"

def plot_fig04():
    t_days = np.linspace(0, 730, 1000)   # 0–24 months
    sol = phase_evolution(t_days)

    fig, ax = plt.subplots(figsize=(5.5, 3.5))

    ax.stackplot(
        t_days / 30.4375,
        sol["x_HAp"],
        sol["x_OCP"],
        sol["x_ACP"],
        labels=["HAp (slow dissolving)", "OCP (moderate)", "ACP (fast dissolving)"],
        colors=["#d62728", "#ff7f0e", "#1f77b4"],
        alpha=0.85,
    )

    # 6-month marker
    ax.axvline(6.0, color="black", lw=1.5, ls="--")
    ax.text(6.2, 0.92, "6 months", fontsize=8, va="top")

    # Phase fractions at 6 months
    idx6 = int(np.argmin(np.abs(t_days - 180)))
    x_acp6 = sol["x_ACP"][idx6]
    x_ocp6 = sol["x_OCP"][idx6]
    x_hap6 = sol["x_HAp"][idx6]
    ax.text(6.3, 0.03,
            f"@6 mo: ACP {x_acp6:.0%} | OCP {x_ocp6:.0%} | HAp {x_hap6:.0%}",
            fontsize=7, color="white", fontweight="bold")

    ax.set_xlabel("Storage duration (months)", fontsize=9)
    ax.set_ylabel("Phase fraction", fontsize=9)
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 3, 6, 9, 12, 18, 24])
    ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9)
    ax.set_title(
        "ACP → OCP → HAp transformation during storage\n"
        "T = −20 °C, pH 7.81 (0% CO₂ loss), glycerol 15%",
        fontsize=8.5,
    )

    fig.tight_layout()
    out = FIGURES_DIR / "fig04_phase_evolution.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  Saved: {out}")

if __name__ == "__main__":
    plot_fig04()
