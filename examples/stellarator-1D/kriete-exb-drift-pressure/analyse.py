#!/usr/bin/env python3
"""Analyse runs of the binormal ExB drift case with evolving pressure.

Usage:
    python3 analyse.py RUN_DIR [RUN_DIR ...] [--plot out.png] [--labels a,b,...]

For each run prints:
  - steadiness (max relative change of n, Te, Ti over the last output)
  - upstream and target temperatures, target densities and their ratio
  - u / c_s at both targets, where u = v_par + v_E/Theta (generalised Bohm)
  - particle balance: total source vs flux through the targets
  - energy balance: (sheath losses + kinetic energy given to source particles)
    / (power input + ExB work)
  - drift strength gamma = v_E / (2 c_s Theta), with c_s at the mid-point

Requires a single-processor run with diagnose = true for the sheath and drift
components (as in BOUT.inp).
"""

import argparse
import glob
import os

import numpy as np
from netCDF4 import Dataset

QE = 1.602176634e-19
MP = 1.67262192595e-27


def conv(ds, name, it=-1):
    var = ds.variables[name]
    data = np.array(var[it] if "t" in var.dimensions else var[:], dtype=float)
    return data * float(getattr(var, "conversion", 1.0))


def load(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "BOUT.dmp.*.nc")))
    if len(files) != 1:
        raise RuntimeError(f"{run_dir}: expected one dump file, found {len(files)}")
    r = {}
    with Dataset(files[0]) as ds:
        myg = int(ds.variables["MYG"][:])
        ny = ds.variables["Ni"].shape[2] - 2 * myg
        sl = slice(myg, myg + ny)

        def prof(name, it=-1):
            return conv(ds, name, it)[0, sl, 0]

        rho_s0 = float(ds.variables["rho_s0"][:])
        dy = np.array(ds.variables["dy"][:])[0, sl]
        g_22 = np.array(ds.variables["g_22"][:])[0, sl] * rho_s0**2
        J = np.array(ds.variables["J"][:])[0, sl]
        # Parallel cell length [m] and cell volume per unit flux-tube area [m]
        ds_par = dy * np.sqrt(g_22)
        r["ds"] = ds_par
        r["s"] = np.cumsum(ds_par) - 0.5 * ds_par
        r["length"] = np.sum(ds_par)
        # Volume element per unit cross-section (J / sqrt(g_22) is the area)
        r["dV"] = ds_par * (J / J.mean())

        r["n"] = prof("Ni")
        r["Te"] = prof("Te")
        r["Ti"] = prof("Ti")
        r["v"] = prof("NVi") / (MP * r["n"])  # AA = 1
        r["w"] = conv(ds, "w_binormal_exb_drift")[0, sl, 0]
        r["work_par"] = prof("Ei_VgradP") + prof("Ee_VgradP")  # v_par dp/ds [W/m^3]

        changes = []
        for name in ["Ni", "Pe", "Pi"]:
            now = prof(name)
            before = prof(name, -2)
            changes.append(np.max(np.abs(now - before) / np.abs(now)))
        r["rel_change"] = max(changes)
        r["nt"] = ds.variables["Ni"].shape[0]

        # Sources [m^-3 s^-1] and [W m^-3]
        r["Sn"] = prof("Si_src")
        r["Spower"] = 1.5 * (prof("Pi_src") + prof("Pe_src"))  # Pa/s -> W/m^3
        # Sheath sinks (in the cells next to the targets)
        r["Esheath"] = prof("Ei_sheath") + prof("Ee_sheath")  # W/m^3, negative
        r["Ssheath"] = prof("Si_sheath")  # m^-3 s^-1 (sign differs by target)
    return r


def analyse(r):
    cs = np.sqrt((r["Ti"] + r["Te"]) * QE / MP)
    mid = len(cs) // 2
    res = {}
    res["gamma"] = float(np.mean(r["w"])) / (2 * cs[mid])
    res["u_lower"] = (r["v"][0] + r["w"][0]) / cs[0]
    res["u_upper"] = (r["v"][-1] + r["w"][-1]) / cs[-1]
    res["ratio"] = r["n"][-1] / r["n"][0]
    res["Te_up"] = np.max(r["Te"])
    res["Te_lower"] = r["Te"][0]
    res["Te_upper"] = r["Te"][-1]
    res["Ti_lower"] = r["Ti"][0]
    res["Ti_upper"] = r["Ti"][-1]

    # Particle balance: source vs. target outflow n |u| (sheath diagnostic)
    particles_in = np.sum(r["Sn"] * r["dV"])
    particles_out = np.sum(np.abs(r["Ssheath"]) * r["dV"])
    res["particle_balance"] = particles_out / particles_in

    # Energy balance. Energy enters through the power source and the ExB
    # work term w dp/ds (a volumetric exchange; for constant w its integral is
    # w (p_upper - p_lower), using the cells next to the targets as an
    # approximation of the face values). It leaves through the sheaths
    # (total gamma T n u) and as the kinetic energy given to source particles,
    # which are injected at rest: Hermes-3 does not return this to heat.
    power_in = np.sum(r["Spower"] * r["dV"])
    power_out = -np.sum(r["Esheath"] * r["dV"])
    p_total = r["n"] * (r["Te"] + r["Ti"]) * QE
    work_exb = float(np.mean(r["w"])) * (p_total[-1] - p_total[0])
    dilution = np.sum(0.5 * MP * r["v"] ** 2 * r["Sn"] * r["dV"])
    res["power_in"] = power_in
    res["work_exb"] = work_exb
    res["dilution"] = dilution
    res["sheath_over_in"] = power_out / power_in
    res["energy_balance"] = (power_out + dilution) / (power_in + work_exb)

    # Kriete's isothermal prediction for the target density ratio with the
    # same drift strength
    res["ratio_kriete"] = (1 + res["gamma"]) / (1 - res["gamma"])
    res["rel_change"] = r["rel_change"]
    return res


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--plot", default=None)
    parser.add_argument("--labels", default=None)
    args = parser.parse_args()

    labels = args.labels.split(",") if args.labels else args.runs
    loaded = []
    for run_dir, label in zip(args.runs, labels):
        r = load(run_dir)
        res = analyse(r)
        loaded.append((label, r, res))
        print(
            f"{label}: gamma={res['gamma']:+.3f} "
            f"| Te upstream {res['Te_up']:.1f} eV, targets {res['Te_lower']:.1f} / "
            f"{res['Te_upper']:.1f} eV | Ti targets {res['Ti_lower']:.1f} / "
            f"{res['Ti_upper']:.1f} eV "
            f"| n_up/n_lo={res['ratio']:.3f} (Kriete isothermal {res['ratio_kriete']:.3f}) "
            f"| u/c_s {res['u_lower']:+.3f} {res['u_upper']:+.3f} "
            f"| particles out/in {res['particle_balance']:.4f} "
            f"| sheath/input {res['sheath_over_in']:.4f}, "
            f"ExB work {100 * res['work_exb'] / res['power_in']:+.1f}%, "
            f"dilution {100 * res['dilution'] / res['power_in']:.1f}% of input "
            f"-> energy balance {res['energy_balance']:.4f} "
            f"| last rel. change {res['rel_change']:.1e} ({r['nt']} outputs)"
        )

    if args.plot:
        plot(loaded, args.plot)


SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
INK = "#3b3b38"
MUTED = "#8a8a84"


def plot(loaded, filename):
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
    for k, (label, r, _) in enumerate(loaded):
        c = SERIES_COLORS[k % len(SERIES_COLORS)]
        x = r["s"] / r["length"] * 2 - 1
        cs = np.sqrt((r["Ti"] + r["Te"]) * QE / MP)
        axes[0].plot(x, r["n"] / 1e19, color=c, lw=2, label=label)
        axes[1].plot(x, r["Te"], color=c, lw=2, label=f"{label}: Te")
        axes[1].plot(x, r["Ti"], color=c, lw=1.2, ls="--", label=f"{label}: Ti")
        axes[2].plot(x, r["v"] / cs, color=c, lw=2, label=f"{label}: v∥/c_s")
        axes[2].plot(
            x,
            (r["v"] + r["w"]) / cs,
            color=c,
            lw=1.2,
            ls="--",
            label=f"{label}: (v∥ + v_E/Θ)/c_s",
        )
    axes[0].set_ylabel("Density [10¹⁹ m⁻³]")
    axes[1].set_ylabel("Temperature [eV]  (solid Te, dashed Ti)")
    axes[2].set_ylabel("Mach number  (solid v∥, dashed v∥ + v_E/Θ)")
    for y in (-1, 0, 1):
        axes[2].axhline(y, color=MUTED, lw=0.6)
    for ax in axes:
        ax.set_xlabel("y / y_target   (lower target ← → upper target)")
        ax.grid(True, color="#e6e6e3", lw=0.6)
        ax.set_axisbelow(True)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle(
        "1D island SOL with binormal E×B drift and evolving pressure (Hermes-3)",
        color=INK,
    )
    fig.savefig(filename, dpi=150)
    print(f"Saved {filename}")


if __name__ == "__main__":
    main()
