#!/usr/bin/env python3
"""Plot target conditions against the recycling fraction R for scans at
different drift velocities.

Usage:
    python3 plot_scan.py OUT.png RUN_DIR [RUN_DIR ...]

Each run directory must contain a single-processor run of BOUT.inp in this
directory (any R and drift:vE). Runs are grouped by vE and sorted by R.
"""

import sys

import numpy as np
from analyse import analyse, load, read_option

SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK = "#3b3b38"
MUTED = "#8a8a84"


def main():
    out, runs = sys.argv[1], sys.argv[2:]
    groups = {}
    for run_dir in runs:
        vE = float(read_option(run_dir, "drift", "vE"))
        r = load(run_dir)
        res = analyse(r)
        groups.setdefault(vE, []).append((r["R"], res))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 9,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for k, vE in enumerate(sorted(groups)):
        c = SERIES_COLORS[k % len(SERIES_COLORS)]
        data = sorted(groups[vE], key=lambda x: x[0])
        pumped = np.array([1 - R for R, _ in data])

        def get(key, data=data):
            return np.array([res[key] for _, res in data])

        label = f"v_E = {vE:g} m/s"
        for ax, lower, upper, scale in [
            (axes[0], "Te_lower", "Te_upper", 1.0),
            (axes[1], "flux_lower", "flux_upper", 1e-24),
        ]:
            ax.plot(
                pumped,
                get(upper) * scale,
                color=c,
                lw=2,
                marker="o",
                ms=4,
                label=f"{label}: upper target",
            )
            ax.plot(
                pumped,
                get(lower) * scale,
                color=c,
                lw=1.2,
                ls="--",
                marker="o",
                ms=4,
                mfc="white",
                label=f"{label}: lower target",
            )
        axes[2].plot(pumped, get("ratio"), color=c, lw=2, marker="o", ms=4, label=label)
    axes[0].set_ylabel("Target electron temperature [eV]")
    axes[1].set_ylabel("Flux into target [10²⁴ m⁻² s⁻¹]")
    axes[2].set_ylabel("Target density ratio  n_upper / n_lower")
    axes[2].set_yscale("log")
    for ax in axes:
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("Pumped fraction 1 − R   (higher recycling →)")
        ax.grid(True, color="#e6e6e3", lw=0.6)
        ax.set_axisbelow(True)
    axes[0].legend(fontsize=7.5, frameon=False)
    axes[2].legend(fontsize=7.5, frameon=False)
    fig.suptitle(
        "Recycling scan at fixed core source and power: the drift (towards the "
        "upper target) makes the targets increasingly asymmetric",
        color=INK,
    )
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
