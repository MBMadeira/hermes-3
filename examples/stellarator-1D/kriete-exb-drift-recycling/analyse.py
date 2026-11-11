#!/usr/bin/env python3
"""Analyse runs of the binormal ExB drift case with neutrals and recycling.

Usage:
    python3 analyse.py RUN_DIR [RUN_DIR ...] [--plot out.png] [--labels a,b,...]

For each run prints:
  - steadiness (max relative change of the evolved fields over the last output)
  - upstream and target temperatures and densities, and the target density ratio
  - u / c_Bohm at both targets, where u = v_par + v_E/Theta (generalised Bohm)
  - flux into each target, n u, and its ratio
  - particle balance: core source + ionisation - recombination vs. target flux,
    and core source vs. pumped flux (1 - R) * target flux
  - energy balance: (sheath losses + radiation + kinetic energy given to core
    source particles) / (power input + ExB work + energy of recycled atoms)

Requires a single-processor run with diagnose = true for the sheath, drift,
recycling and reactions components (as in BOUT.inp).
"""

import argparse
import glob
import os

import numpy as np
from netCDF4 import Dataset

QE = 1.602176634e-19
MP = 1.67262192595e-27

ION = "d+"
NEUTRAL = "d"


def conv(ds, name, it=-1):
    var = ds.variables[name]
    data = np.array(var[it] if "t" in var.dimensions else var[:], dtype=float)
    return data * float(getattr(var, "conversion", 1.0))


def read_option(run_dir, section, key):
    """Read an evaluated option from BOUT.settings"""
    current = ""
    with open(os.path.join(run_dir, "BOUT.settings")) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if line.startswith("["):
                current = line.strip("[]")
            elif "=" in line and current == section:
                name, value = (s.strip() for s in line.split("=", 1))
                if name == key:
                    return value
    raise KeyError(f"{section}:{key}")


def load(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "BOUT.dmp.*.nc")))
    if len(files) != 1:
        raise RuntimeError(f"{run_dir}: expected one dump file, found {len(files)}")
    r = {}
    with Dataset(files[0]) as ds:
        myg = int(ds.variables["MYG"][:])
        ny = ds.variables[f"N{ION}"].shape[2] - 2 * myg
        sl = slice(myg, myg + ny)

        def prof(name, it=-1):
            return conv(ds, name, it)[0, sl, 0]

        rho_s0 = float(ds.variables["rho_s0"][:])
        dy = np.array(ds.variables["dy"][:])[0, sl]
        g_22 = np.array(ds.variables["g_22"][:])[0, sl] * rho_s0**2
        J = np.array(ds.variables["J"][:])[0, sl]
        ds_par = dy * np.sqrt(g_22)  # Parallel cell length [m]
        r["s"] = np.cumsum(ds_par) - 0.5 * ds_par
        r["length"] = np.sum(ds_par)
        r["dV"] = ds_par * (J / J.mean())  # Volume per unit flux-tube area [m]
        r["AA"] = float(read_option(run_dir, ION, "AA"))
        r["R"] = float(read_option(run_dir, ION, "target_recycle_multiplier"))

        r["n"] = prof(f"N{ION}")
        r["Te"] = prof("Te")
        r["Ti"] = prof(f"T{ION}")
        r["nn"] = prof(f"N{NEUTRAL}")
        r["Tn"] = prof(f"T{NEUTRAL}")
        r["v"] = prof(f"V{ION}")
        r["w"] = conv(ds, "w_binormal_exb_drift")[0, sl, 0]

        changes = []
        for name in [f"N{ION}", f"P{ION}", "Pe", f"N{NEUTRAL}", f"P{NEUTRAL}"]:
            now = prof(name)
            before = prof(name, -2)
            changes.append(np.max(np.abs(now - before) / np.abs(now)))
        r["rel_change"] = max(changes)
        r["nt"] = ds.variables[f"N{ION}"].shape[0]

        r["Sn"] = prof(f"S{ION}_src")  # Core particle source [m^-3 s^-1]
        r["Spower"] = 1.5 * (prof(f"P{ION}_src") + prof("Pe_src"))  # [W m^-3]
        r["Ssheath"] = prof(f"S{ION}_sheath")  # Sign differs by target
        r["Esheath"] = prof(f"E{ION}_sheath") + prof("Ee_sheath")  # Negative
        r["Siz"] = prof(f"S{ION}_iz")
        r["Srec"] = prof(f"S{ION}_rec")  # Negative (ion sink)
        r["Srecycle"] = prof(f"S{NEUTRAL}_target_recycle")
        r["Erecycle"] = prof(f"E{NEUTRAL}_target_recycle")
        r["Rad"] = prof(f"R{ION}_ex") + prof(f"R{ION}_rec")  # Negative (loss)
        # Other radiation and reaction losses, if any: fixed fraction
        # impurities, and electron energy losses of other reactions (e.g.
        # impurity ionisation and recombination)
        for name in ds.variables:
            if name.startswith("Rfixed_fraction") or (
                name.startswith("R")
                and name.endswith(("_ex", "_rec"))
                and name not in (f"R{ION}_ex", f"R{ION}_rec")
                and "t" in ds.variables[name].dimensions
            ):
                r["Rad"] = r["Rad"] + prof(name)
    return r


def analyse(r):
    mi = r["AA"] * MP
    cs = np.sqrt((r["Ti"] + r["Te"]) * QE / mi)
    mid = len(cs) // 2
    dV = r["dV"]
    res = {}
    res["gamma"] = float(np.mean(r["w"])) / (2 * cs[mid])
    res["u_lower"] = (r["v"][0] + r["w"][0]) / cs[0]
    res["u_upper"] = (r["v"][-1] + r["w"][-1]) / cs[-1]
    res["ratio"] = r["n"][-1] / r["n"][0]
    for key in ["Te", "Ti"]:
        res[f"{key}_up"] = r[key][mid]
        res[f"{key}_lower"] = r[key][0]
        res[f"{key}_upper"] = r[key][-1]
    res["n_up"] = r["n"][mid]
    res["n_lower"] = r["n"][0]
    res["n_upper"] = r["n"][-1]
    res["nn_lower"] = r["nn"][0]
    res["nn_upper"] = r["nn"][-1]

    # Target fluxes from the sheath diagnostic, in the cells next to each target
    res["flux_lower"] = abs(r["Ssheath"][0]) * dV[0]
    res["flux_upper"] = abs(r["Ssheath"][-1]) * dV[-1]
    flux_targets = res["flux_lower"] + res["flux_upper"]
    res["flux_targets"] = flux_targets

    core = np.sum(r["Sn"] * dV)
    ionisation = np.sum(r["Siz"] * dV)
    recombination = -np.sum(r["Srec"] * dV)
    recycled = np.sum(r["Srecycle"] * dV)
    res["core"] = core
    res["ionisation"] = ionisation
    # Ions: core + ionisation - recombination = flux into the targets
    res["ion_balance"] = flux_targets / (core + ionisation - recombination)
    # Recycling: recycled flux = R * flux into the targets
    res["recycled_over_R_flux"] = recycled / (r["R"] * flux_targets)
    # Total: core source = pumped flux (1 - R) * flux into the targets
    res["pump_balance"] = (1 - r["R"]) * flux_targets / core

    # Energy balance (all species, internal + kinetic energy). The reactions
    # conserve energy apart from the electron losses (radiation and ionisation
    # potential, R diagnostics). The neutrals have no-flow boundaries.
    power_in = np.sum(r["Spower"] * dV)
    sheath = -np.sum(r["Esheath"] * dV)
    radiation = -np.sum(r["Rad"] * dV)
    recycle_energy = np.sum(r["Erecycle"] * dV)
    p_total = r["n"] * (r["Te"] + r["Ti"]) * QE
    work_exb = float(np.mean(r["w"])) * (p_total[-1] - p_total[0])
    dilution = np.sum(0.5 * mi * r["v"] ** 2 * r["Sn"] * dV)
    res["power_in"] = power_in
    res["sheath"] = sheath
    res["radiation"] = radiation
    res["work_exb"] = work_exb
    res["recycle_energy"] = recycle_energy
    res["dilution"] = dilution
    res["energy_balance"] = (sheath + radiation + dilution) / (
        power_in + work_exb + recycle_energy
    )
    res["ratio_kriete"] = (1 + res["gamma"]) / (1 - res["gamma"])
    res["rel_change"] = r["rel_change"]
    return res


def summary(label, r, res):
    return (
        f"{label}: gamma={res['gamma']:+.3f} | upstream n={res['n_up']:.3g} m^-3, "
        f"Te={res['Te_up']:.1f} eV, Ti={res['Ti_up']:.1f} eV\n"
        f"    targets (lower / upper): n {res['n_lower']:.3g} / {res['n_upper']:.3g} "
        f"(ratio {res['ratio']:.3f}, Kriete isothermal {res['ratio_kriete']:.3f}), "
        f"Te {res['Te_lower']:.1f} / {res['Te_upper']:.1f} eV, "
        f"Ti {res['Ti_lower']:.1f} / {res['Ti_upper']:.1f} eV, "
        f"n_atom {res['nn_lower']:.3g} / {res['nn_upper']:.3g}\n"
        f"    flux into targets {res['flux_lower']:.3g} / {res['flux_upper']:.3g} "
        f"m^-2 s^-1 (upper/lower {res['flux_upper'] / res['flux_lower']:.3f}), "
        f"u/c_s {res['u_lower']:+.3f} {res['u_upper']:+.3f}\n"
        f"    particles: target flux / (core + iz - rec) {res['ion_balance']:.4f}, "
        f"recycled / (R flux) {res['recycled_over_R_flux']:.4f}, "
        f"(1-R) flux / core {res['pump_balance']:.4f}\n"
        f"    energy: input {res['power_in']:.3g} W/m^2, ExB work "
        f"{100 * res['work_exb'] / res['power_in']:+.1f}%, sheaths "
        f"{100 * res['sheath'] / res['power_in']:.1f}%, radiation + ionisation "
        f"{100 * res['radiation'] / res['power_in']:.1f}% of input "
        f"-> balance {res['energy_balance']:.4f}\n"
        f"    last rel. change {res['rel_change']:.1e} ({r['nt']} outputs)"
    )


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
        print(summary(label, r, res))

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
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    axes = axes.flatten()
    for k, (label, r, _) in enumerate(loaded):
        c = SERIES_COLORS[k % len(SERIES_COLORS)]
        x = r["s"] / r["length"] * 2 - 1
        cs = np.sqrt((r["Ti"] + r["Te"]) * QE / (r["AA"] * MP))
        axes[0].plot(x, r["n"] / 1e19, color=c, lw=2, label=label)
        axes[1].plot(x, r["Te"], color=c, lw=2, label=f"{label}: Te")
        axes[1].plot(x, r["Ti"], color=c, lw=1.2, ls="--", label=f"{label}: Ti")
        axes[2].semilogy(x, r["nn"], color=c, lw=2, label=label)
        axes[3].plot(x, r["v"] / cs, color=c, lw=2, label=f"{label}: v∥/c_s")
        axes[3].plot(
            x,
            (r["v"] + r["w"]) / cs,
            color=c,
            lw=1.2,
            ls="--",
            label=f"{label}: (v∥ + v_E/Θ)/c_s",
        )
    axes[0].set_ylabel("Ion density [10¹⁹ m⁻³]")
    axes[1].set_ylabel("Temperature [eV]  (solid Te, dashed Ti)")
    axes[2].set_ylabel("Atom density [m⁻³]")
    axes[3].set_ylabel("Mach number  (solid v∥, dashed v∥ + v_E/Θ)")
    for y in (-1, 0, 1):
        axes[3].axhline(y, color=MUTED, lw=0.6)
    for ax in axes:
        ax.set_xlabel("y / y_target   (lower target ← → upper target)")
        ax.grid(True, color="#e6e6e3", lw=0.6)
        ax.set_axisbelow(True)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle(
        "1D island SOL with binormal E×B drift, neutrals and recycling (Hermes-3)",
        color=INK,
    )
    fig.savefig(filename, dpi=150)
    print(f"Saved {filename}")


if __name__ == "__main__":
    main()
