#!/usr/bin/env python3
"""Compare Hermes-3 runs of the binormal ExB drift benchmark with Kriete's
analytic solution.

Usage:
    python3 compare.py RUN_DIR [RUN_DIR ...] [--plot out.png]

Each RUN_DIR must contain BOUT.dmp.*.nc from a single-processor run of
BOUT.inp in this directory (any gamma, Theta, resolution).
"""

import argparse
import glob
import os

import analytic
import numpy as np
from netCDF4 import Dataset

QE = 1.602176634e-19
MP = 1.67262192595e-27


def conv(ds, name):
    """Variable converted to SI units using its 'conversion' attribute."""
    var = ds.variables[name]
    data = np.array(var[:], dtype=float)
    return data * float(getattr(var, "conversion", 1.0))


def load(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "BOUT.dmp.*.nc")))
    if len(files) != 1:
        raise RuntimeError(f"{run_dir}: expected one dump file, found {len(files)}")
    with Dataset(files[0]) as ds:
        myg = int(ds.variables["MYG"][:])
        Ni = conv(ds, "Ni")[-1, 0, :, 0]  # Last time, x=0, z=0 [m^-3]
        NVi = conv(ds, "NVi")[-1, 0, :, 0]  # [kg m^-2 s^-1]
        Ni_prev = conv(ds, "Ni")[-2, 0, :, 0]
        rho_s0 = float(ds.variables["rho_s0"][:])
        # dy is not normalised; g_22 is normalised by rho_s0^2, so the
        # parallel length of a cell is dy * sqrt(g_22) * rho_s0 [m]
        dy = np.array(ds.variables["dy"][:])[0, :]
        g_22 = np.array(ds.variables["g_22"][:])[0, :] * rho_s0**2
        Ti = float(conv(ds, "Ti").ravel()[0])  # [eV]
        Te = float(conv(ds, "Te").ravel()[0])
        w = conv(ds, "w_binormal_exb_drift")[0, :, 0]  # v_E / Theta [m/s]
        AA = 1.0

    ny = Ni.size - 2 * myg
    interior = slice(myg, myg + ny)
    # Parallel length of each cell and its centre position
    ds_par = (dy * np.sqrt(g_22))[interior]
    s_face = np.concatenate([[0.0], np.cumsum(ds_par)])
    s = 0.5 * (s_face[1:] + s_face[:-1])
    length = s_face[-1]

    n = Ni[interior]
    v = NVi[interior] / (AA * MP * n)
    cs = np.sqrt((Ti + Te) * QE / (AA * MP))

    # Values in the cells next to the targets. Note: guard cells in the dump
    # file hold the solver's values, not those set by the sheath component.
    rel_change = np.max(np.abs(Ni - Ni_prev)[interior] / Ni[interior])

    settings = read_settings(os.path.join(run_dir, "BOUT.settings"))

    return {
        "s": s,
        "length": length,
        "n": n,
        "v": v,
        "cs": cs,
        "w": w[interior],
        "n_lower": n[0],
        "n_upper": n[-1],
        "v_lower": v[0],
        "v_upper": v[-1],
        "rel_change": rel_change,
        "ny": ny,
        "xpoint_frac": float(settings["mesh"]["y_xpoint"])
        / float(settings["mesh"]["y_target"]),
        "nmax_input": float(settings["kriete"]["nmax"]),
    }


def read_settings(path):
    """Minimal parser for BOUT.settings: {section: {key: value_string}}."""
    sections = {"": {}}
    current = ""
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                current = line[1:-1]
                sections.setdefault(current, {})
            elif "=" in line:
                key, value = line.split("=", 1)
                sections[current][key.strip()] = value.strip()
    return sections


def zero_crossing(x, f):
    """Linear interpolation of the first sign change of f(x)."""
    idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    return x[i] - f[i] * (x[i + 1] - x[i]) / (f[i + 1] - f[i])


def analyse(run):
    """Compare one run with the analytic solution."""
    cs = run["cs"]
    w = run["w"]
    gamma = float(np.mean(w)) / (2 * cs)  # gamma = v_E / (2 c_s Theta) = w / (2 c_s)

    half = run["length"] / 2
    s_x = run["xpoint_frac"] * half
    ystar = (run["s"] - half) / s_x

    # Analytic peak density, set in the input through the source:
    # n_max = 2 S s_x / c_s
    n_max = run["nmax_input"]

    n_an, M_an = analytic.solution(ystar, gamma)
    nstar = run["n"] / n_max
    M = run["v"] / cs

    res = {
        "gamma": gamma,
        "ystar": ystar,
        "nstar": nstar,
        "M": M,
        "n_an": n_an,
        "M_an": M_an,
    }
    res["err_n_max"] = np.max(np.abs(nstar - n_an))
    res["err_M_max"] = np.max(np.abs(M - M_an))
    res["ratio"] = run["n_upper"] / run["n_lower"]
    res["ratio_an"] = analytic.target_density_ratio(gamma)
    res["M_upper"] = run["v_upper"] / cs
    res["M_lower"] = run["v_lower"] / cs
    res["u_upper"] = (run["v_upper"] + np.mean(w)) / cs
    res["u_lower"] = (run["v_lower"] + np.mean(w)) / cs
    res["y_stag"] = zero_crossing(ystar, M)
    res["y_stag_an"] = analytic.stagnation_point(gamma)
    res["y_vy0"] = zero_crossing(ystar, M + 2 * gamma)
    res["y_vy0_an"] = analytic.zero_poloidal_flux_point(gamma)
    # Particle balance: total source / c_s = n_max, total target outflow
    # (n |u| at both targets) / c_s should equal it in steady state
    res["outflow_over_nmax"] = (
        run["n_upper"] * abs(res["u_upper"]) + run["n_lower"] * abs(res["u_lower"])
    ) / n_max
    res["xpoint_frac"] = run["xpoint_frac"]
    res["rel_change"] = run["rel_change"]
    return res


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--plot", default=None)
    args = parser.parse_args()

    results = []
    for run_dir in args.runs:
        run = load(run_dir)
        res = analyse(run)
        res["name"] = run_dir
        res["ny"] = run["ny"]
        results.append(res)
        print(
            f"{run_dir}: ny={run['ny']} gamma={res['gamma']:+.3f} "
            f"| max|dn*|={res['err_n_max']:.4f} max|dM|={res['err_M_max']:.4f} "
            f"| n_up/n_lo={res['ratio']:.4f} (analytic {res['ratio_an']:.4f}) "
            f"| y*_stag={res['y_stag']:+.4f} ({res['y_stag_an']:+.4f}) "
            f"| y*(v_y=0)={res['y_vy0']:+.4f} ({res['y_vy0_an']:+.4f}) "
            f"| u/c_s at targets {res['u_lower']:+.3f} {res['u_upper']:+.3f} "
            f"| outflow/source={res['outflow_over_nmax']:.4f} "
            f"| last-output rel. change={res['rel_change']:.1e}"
        )

    if args.plot:
        plot(results, args.plot)


# Fixed categorical colour order (validated for colour-vision deficiency);
# each series also carries a direct label and a legend entry
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
INK = "#3b3b38"
MUTED = "#8a8a84"


def plot(results, filename):
    """Density and Mach number profiles, Hermes-3 (markers) vs analytic (lines)."""
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
    results = sorted(results, key=lambda r: r["gamma"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    yend = max(1 / r["xpoint_frac"] for r in results)
    yfine = np.linspace(-yend, yend, 1001)

    for k, res in enumerate(results):
        color = SERIES_COLORS[k % len(SERIES_COLORS)]
        label = f"γ = {res['gamma']:+.2f}"
        n_an, M_an = analytic.solution(yfine, res["gamma"])
        # Thin out markers so the analytic line stays visible
        step = max(1, len(res["ystar"]) // 40)
        for ax, sim, an in (
            (axes[0], res["nstar"], n_an),
            (axes[1], res["M"], M_an),
        ):
            ax.plot(yfine, an, color=color, lw=2, label=label)
            ax.plot(
                res["ystar"][::step],
                sim[::step],
                "o",
                color=color,
                ms=4,
                mec="white",
                mew=0.8,
            )
        # Direct label at the right-hand end of the Mach number curve
        axes[1].annotate(
            label,
            (yfine[-1], M_an[-1]),
            xytext=(4, 0),
            textcoords="offset points",
            va="center",
            fontsize=8,
            color=INK,
        )

    for ax in axes:
        for x in (-1, 1):
            ax.axvline(x, color=MUTED, lw=0.8, ls=":")
        ax.set_xlabel("y / y_X-point   (lower target ← → upper target)")
        ax.grid(True, color="#e6e6e3", lw=0.6)
        ax.set_axisbelow(True)
    axes[1].axhline(0, color=MUTED, lw=0.8)
    axes[1].set_xlim(right=yend * 1.18)
    axes[0].set_ylabel("Density  n / n_max")
    axes[1].set_ylabel("Parallel Mach number  v∥ / c_s")

    # Legend: colours identify gamma; line vs marker identifies the source
    handles, labels = axes[0].get_legend_handles_labels()
    handles += [
        plt.Line2D([], [], color=MUTED, lw=2),
        plt.Line2D([], [], color=MUTED, marker="o", ls="", ms=4, mec="white"),
    ]
    labels += ["Kriete et al. analytic", "Hermes-3"]
    axes[0].legend(
        handles, labels, fontsize=8, loc="lower center", frameon=False, ncol=2
    )
    fig.suptitle(
        "Binormal E×B drift in a 1D island SOL: Hermes-3 vs Kriete et al., "
        "Nucl. Fusion 63 (2023) 026022",
        color=INK,
    )
    fig.savefig(filename, dpi=150)
    print(f"Saved {filename}")


if __name__ == "__main__":
    main()
