# Binormal E×B drift with neutrals and recycling

Extends `../kriete-exb-drift-pressure` with deuterium atoms and recycling at
both targets, using the existing Hermes-3 components:

- `recycling`: d+ is recycled as d at both targets. With the option
  `exb_drift = binormal_exb_drift`, the recycled flux is the total flux into
  each target, n (v∥ + v_E/Θ), and fast recycling (not used here) uses
  γ_i T n (v∥ + v_E/Θ). Without the option, the drift's share of the target
  flux would be lost to the pump.
- `reactions`: ionisation, recombination and charge exchange (AMJUEL rates).
- `neutral_parallel_diffusion` with `dneut = 1/Θ²`: cross-field neutral
  diffusion in the island poloidal direction, projected onto the field line
  coordinate s. This is the island analogue of (B/B_pol)² in a tokamak. The
  binormal drift does not act on neutrals.

The core particle source (2×10²³ m⁻² s⁻¹) and power (2 × 6×10⁷ W/m², split
equally between ions and electrons) are fixed and uniform along the flux tube,
as in the pressure example. A fraction 1 − R of the target flux is pumped
(R = `d+:target_recycle_multiplier`, default 0.9), so in steady state the flux
into the two targets is Γ_core / (1 − R). Θ = 10⁻³, L∥ = 1000 m, ny = 200,
packed towards both targets (1 m cells at the targets, 9 m at the mid-point).

## Running

```bash
hermes-3 -d .            # about 1 minute on one core
python3 analyse.py . --plot profiles.png
```

The neutral diffusion along s is ~10⁹ m²/s, which makes the problem very
stiff. The input therefore uses CVODE with the band-block-diagonal
preconditioner (`solver:cvode_precon_method = bbd`, bandwidths 21 = 3 × the
7 evolving variables). Without a preconditioner CVODE fails in the initial
transient; with it the run takes about a minute. The preconditioner
only couples cells on the same processor, so in parallel it is much less
effective: run this case on one processor. Above R = 0.95 without drift, and
with v_E from the sheath potential, CVODE stalls; use the Newton solver
instead (see "Numerics" below; needs BOUT++ built with PETSc).

`analyse.py` prints the target conditions, u/c_s at both targets and the
particle and energy balance:

- ions: flux into the targets = core source + ionisation − recombination
- recycling: recycled flux = R × flux into the targets
- pumping: core source = (1 − R) × flux into the targets
- energy: sheath losses + radiation and ionisation losses + kinetic energy
  given to core source particles = power input + E×B work w (p_upper − p_lower)
  + energy of the recycled atoms (3 eV each). The reactions conserve the total
  (internal + kinetic) energy of the heavy species; the atoms have no-flow
  boundaries.

`plot_scan.py OUT.png RUN_DIR ...` plots the target conditions against 1 − R
for runs at different v_E.

## Results

Recycling scan at fixed core source and power, each run restarted from the
previous one (CVODE: 60 outputs of 5×10⁶/Ω_ci; rows marked N: Newton solver,
converged to a relative change below 10⁻⁶ per output):

| v_E [m/s] | R | n_up [10¹⁹ m⁻³] | T_e lower / upper [eV] | Flux lower / upper [10²⁴ m⁻² s⁻¹] | n_upper / n_lower | u/c_s lower, upper | Energy balance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.90 | 0.86 | 58.8 / 58.8 | 1.00 / 1.00 | 1.00 | −0.95, +0.95 | 1.028 |
| 0 | 0.95 | 1.17 | 28.5 / 28.5 | 2.00 / 2.00 | 1.00 | −0.91, +0.91 | 1.058 |
| 0 | 0.97 (N) | 1.54 | 14.9 / 14.9 | 3.33 / 3.33 | 1.00 | −0.89, +0.89 | 1.057 |
| 0 | 0.98 (N) | 1.86 | 9.6 / 9.6 | 5.00 / 5.00 | 1.00 | −0.87, +0.87 | 1.068 |
| 0 | 0.99 (N) | 2.46 | 3.8 / 3.8 | 10.0 / 10.0 | 1.00 | −0.86, +0.86 | 1.049 |
| 6.6 | 0.90 | 0.86 | 63.7 / 54.7 | 0.90 / 1.10 | 1.34 | −0.95, +0.94 | 1.027 |
| 6.6 | 0.95 | 1.16 | 32.9 / 25.2 | 1.72 / 2.28 | 1.48 | −0.91, +0.91 | 1.055 |
| 6.6 | 0.97 (N) | 1.50 | 18.4 / 13.0 | 2.74 / 3.93 | 1.66 | −0.89, +0.89 | 1.059 |
| 6.6 | 0.98 (N) | 1.81 | 12.2 / 7.8 | 3.83 / 6.17 | 1.98 | −0.88, +0.88 | 1.059 |
| 66 | 0.90 | 0.59 | 126.5 / 44.5 | 0.28 / 1.72 | 11.0 | −0.98, +0.98 | 0.996 |
| 66 | 0.95 | 0.61 | 122.9 / 17.2 | 0.29 / 3.71 | 20.8 | −0.98, +1.64 | 0.993 |
| 66 | 0.97 | 0.62 | 120.0 / 7.9 | 0.30 / 6.37 | 33.3 | −0.98, +2.42 | 0.989 |
| 66 | 0.98 | 0.61 | 120.4 / 3.9 | 0.29 / 9.71 | 54.1 | −0.98, +3.50 | 0.990 |

n_up is the density at the mid-point. Particle balance closes to 10⁻⁵ in all
cases. With v_E = 6.6 and 66 m/s, γ = v_E/(2 c_s Θ) ≈ 0.025 and 0.25 at the
mid-point, for which Kriete's isothermal model (no neutrals) gives
n_upper/n_lower = 1.05 and 1.67.

- Recycling strongly amplifies the drift asymmetry. The drift carries more
  plasma to the upper target, which recycles more, cools and becomes denser,
  which in turn increases the ionisation there. Even v_E = 6.6 m/s makes the
  target fluxes differ by 20–30%. With v_E = 66 m/s nearly all the extra
  recycling flux goes to the upper target: the lower target stays at
  ~120 eV and ~3×10²³ m⁻² s⁻¹ throughout the scan, while the upper target
  cools from 45 to 4 eV.
- At the cold target the drift becomes faster than sound: w = v_E/Θ =
  66 km/s, while c_s ≈ 20 km/s at 4 eV. The drift alone then carries the
  plasma into the target (u > c_s with v∥ ≈ 0), and the Bohm condition no
  longer constrains v∥ there. With constant v_E this is where the model is
  least realistic: in W7-X, v_E falls as the target plasma cools.
- At R = 0.99 (target flux Γ_core/(1 − R) = 2×10²⁵ m⁻² s⁻¹) there is a
  steady state without drift, at 3.8 eV at both targets, with radiation and
  ionisation taking 52% of the input power. With any drift no steady state was
  found: with v_E = 66 m/s the plasma cools to ~1 eV along the whole flux
  tube and the atoms accumulate (CVODE); with v_E = 6.6 m/s, or v_E from the
  sheath potential, the upper target cools below 2 eV and the Newton solver
  makes no further progress.
- The energy balance closes to ~1% with the strong drift and to 3–6% in the
  near-symmetric cases. The residual is a discretisation error at the steep
  target gradients: at R = 0.9, v_E = 0 it falls from 2.8% to 1.6% when ny is
  doubled to 400.
- Without drift, or with the weak drift, CVODE cannot go beyond R = 0.95: it
  stalls at R = 0.97 (from a restart and from scratch, at ny = 100 and 200),
  when the atom density at the mid-point falls below ~2×10¹³ m⁻³. Neither a
  larger Krylov space, tighter absolute tolerances for the atoms, nor atoms at
  fixed zero velocity removed the stall. The Newton solver reaches these steady
  states (rows marked N). CVODE restarted from the converged Newton solution at
  R = 0.97 fails on its first step, so the problem is CVODE's linear solve
  (the banded difference-quotient preconditioner), not the absence of a steady
  state. With the strong drift the hot, tenuous side lets atoms reach the
  mid-point (~10¹⁵ m⁻³) and CVODE works.

Other checks (R = 0.9, v_E = 66 m/s): v_E → −v_E mirrors the solution
(2×10⁻⁶ at ny = 100); ny = 400 gives n_upper/n_lower = 11.07 and a flux ratio
of 6.01 (11.00 and 6.05 at ny = 200); the right-hand side is identical on 1 and
2 processors.

### v_E from the sheath potential

With `binormal_exb_drift:vE_mode = sheath` (see the commented lines in
`BOUT.inp`), v_E = φ/(B λ_φ) follows the target sheath potentials. Here
λ_φ = 1.439 m is chosen so that the R = 0.9 steady state with constant
v_E = 66 m/s is a fixed point (φ = 348 V at the lower target and 127 V at the
upper; mean 238 V). The same R scan (Newton solver from R = 0.95 on):

| v_E mode | R | v_E [m/s] | φ lower / upper [V] | T_e lower / upper [eV] | Flux lower / upper [10²⁴ m⁻² s⁻¹] | n_upper / n_lower | u/c_s upper | Energy balance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| constant | 0.90 | 66 | | 126.5 / 44.5 | 0.28 / 1.72 | 11.0 | +0.98 | 0.996 |
| sheath | 0.90 | 66.0 | 348 / 127 | 126.5 / 44.5 | 0.28 / 1.72 | 11.0 | +0.98 | 0.996 |
| constant | 0.95 | 66 | | 122.9 / 17.2 | 0.29 / 3.71 | 20.8 | +1.64 | 0.993 |
| sheath | 0.95 | 36.5 | 207 / 55 | 72.2 / 19.4 | 0.70 / 3.30 | 8.45 | +0.95 | 1.010 |
| constant | 0.97 | 66 | | 120.0 / 7.9 | 0.30 / 6.37 | 33.3 | +2.42 | 0.989 |
| sheath | 0.97 | 14.0 | 69 / 33 | 24.1 / 11.5 | 2.13 / 4.54 | 2.89 | +0.90 | 1.051 |
| constant | 0.98 | 66 | | 120.4 / 3.9 | 0.29 / 9.71 | 54.1 | +3.50 | 0.990 |
| sheath | 0.98 | 8.0 | 36 / 21 | 13.0 / 7.5 | 3.62 / 6.38 | 2.26 | +0.88 | 1.058 |

- With a constant v_E, raising R cools only the upper target and the
  asymmetry grows without limit. With v_E from the sheath potential, the
  colder upper target lowers φ and hence v_E, which weakens the drift, so the
  lower target receives more plasma and cools as well. The asymmetry then
  falls with R (11 → 2.3), and the drift never reaches the upper target
  supersonically (u ≈ c_s).
- Which potential sets E_r matters. With `sheath_potential = lower`
  (λ_φ = 2.111 m, again 66 m/s at R = 0.9), the hot lower target keeps v_E
  high: 62 and 53 m/s at R = 0.95 and 0.97 (ratio 19.4 and 23.6, lower
  target 119 and 100 eV), close to the constant-v_E results, before the lower
  target cools at R = 0.98 (v_E = 6.4 m/s, ratio 1.95). With
  `sheath_potential = upper` (λ_φ = 0.768 m), v_E falls fastest at first (29.5
  m/s at R = 0.95, ratio 5.8) and ends at 10.4 m/s at R = 0.98 (ratio 2.9).
- At R = 0.99 no steady state was found (the upper target cools below 2 eV;
  see above).

### Other components with the drift

Added one at a time to this case (R = 0.9, v_E = 66 m/s); each runs to a steady
state, with particle balance to 10⁻⁵, energy balance to ~1–2% and the same
right-hand side on 1 and 2 processors:

| Added | Effect |
| --- | --- |
| `braginskii_ion_viscosity` | n_upper/n_lower 11.0 → 7.4: the lower target, with supersonic v∥, is hot and nearly collisionless (Ti 174 → 122 eV). Consider the flux limiter `eta_limit_alpha`. |
| `braginskii_thermal_force` | None with deuterium only (the electron–ion force cancels at zero current); acts on impurities. Must be listed before `electron_force_balance`. |
| `fixed_fraction_nitrogen` (2%) | Radiation 5% of the input; asymmetry 11 → 12. Mirror symmetric. |
| Trace He⁺ (core source, no atoms) | He⁺ target fluxes 1.7:1 (D 4.6:1); He concentration higher at the hot target. Mirror symmetric. |
| `neutral_boundary` on d | Negligible here: the atoms lose 5 kW/m² to the targets (input 120 MW/m²). |
| `simple_pump` on d (R = 1) | Pumped atoms = core source. |
| `binormal_stpm` (D = 1, χ = 3, ν = 1 m²/s) | Asymmetry 11 → 12.7. Note that it also applies D, χ, ν to the atoms. |
| Density and temperature feedback at the mid-point | See `tests/integrated/1D-binormal-exb-drift-feedback`. |

Not usable, independently of the drift (same at v_E = 0):

- `braginskii_electron_viscosity`: NaN at the targets (guard cells of the
  viscosity coefficient are not set). With that patched, its force is
  10⁴–10⁵ times the Braginskii estimate and it does not conserve energy.
- `snb_conduction`: computes the heat flux through the target faces as well,
  on top of the sheath heat flux (energy balance 0.62). Its bug that rescaled
  `dy` on every call is fixed.
- Evolving impurity atoms (helium): CVODE stalls in the initial transient with
  dneut = 1/Θ². The Newton solver progresses, but its timesteps stay at
  10³–10⁴/Ω_ci, far too short to reach a steady state. In the default `afn`
  mode of `neutral_parallel_diffusion` the helium atoms' only collision is
  ionisation (there is no charge exchange), so their diffusion coefficient is
  very large where the plasma is cold; `diffusion_collisions_mode =
  multispecies` makes the timesteps several times longer, but not enough.

### Processor-boundary fix in `braginskii_conduction`

The 1 vs 2 processor check found a bug in `braginskii_conduction`, unrelated to
the drift: the conduction coefficient κ was only communicated between
processors when the heat flux limiter was used. The charge exchange and
ionisation collision frequencies from `reactions` are only calculated in the
domain interior, so for ions and atoms κ was wrong in the processor guard cells,
and so was the conductive flux across processor boundaries. This is invisible in
symmetric cases (no temperature gradient at the mid-point, where the processor
boundary is on 2 processors) and made the drift cases fail on 2 processors.
κ is now always communicated; single-processor results are unchanged.

## Numerics: the Newton solver

With BOUT++ built with PETSc (`-DBOUT_USE_PETSC=ON`), the backward Euler
solver `beuler` finds the steady states that CVODE cannot reach. Replace the
`[solver]` section of `BOUT.inp` with the commented one below it:

- `equation_form = pseudo_transient`: local pseudo-timesteps, set from the
  local residual, so that the solution moves towards the steady state rather
  than following the time evolution.
- `pc_type = lu`: the preconditioner is the LU factorisation of the
  finite-difference Jacobian, built with colouring. In 1D this is cheap (7 × ny
  unknowns), and it handles the stiff neutral diffusion directly.
- `matrix_free_operator = true` with `ksp_type = gmres`: the Jacobian-vector
  products are evaluated exactly. This is needed with v_E from the sheath
  potential, where v_E couples every cell to the cells next to the targets.
  The coloured Jacobian only knows the local stencil, so with it alone Newton
  converges slowly and the timestep stays near 500/Ω_ci.
- `atol = 1e-11`, `rtol = 1e-8`: with `atol = 1e-7` the solver accepts each
  step without iterating once the residual is below 10⁻⁷, and the solution
  stops about 1% short of the steady state while appearing converged.

Results (one core):

- The CVODE scan steps (60 outputs of 5×10⁶/Ω_ci, about 4 minutes each) take
  8–90 s each and converge to a relative change of 10⁻¹²–10⁻¹⁶ per output.
  From scratch, R = 0.9 with v_E = 66 m/s takes 74 s (CVODE: about a minute)
  and agrees with CVODE to all printed digits.
- R = 0.97–0.99 without drift, R = 0.97 and 0.98 with v_E = 6.6 m/s, and
  v_E from the sheath potential up to R = 0.98, all of which stall CVODE,
  converge.
- On 2 processors, `pc_type = bjacobi` with `sub_pc_type = lu` in a `[petsc]`
  section gives the same steady state as 1 processor to 5×10⁻⁹.
- The full sheath model (`sheath_boundary`, pressure example) at ny = 400,
  which does not settle with CVODE, has a steady state (n_upper/n_lower =
  1.620 at v_E = 66 m/s). CVODE restarted from it stays there (relative change
  10⁻¹¹ per output over 30 outputs).

The test `tests/integrated/1D-binormal-exb-drift-newton` (built only with
PETSc) runs R = 0.97 with v_E from the sheath potential from scratch (about
2 minutes) and checks the steady state, the balances, v_E and a restart on 2
processors. From scratch, the same case converges at ny = 200 (99 s) but made
no progress in 10 minutes at ny = 100: restarting from a nearby solution is the
reliable route.

Where the Newton solver does not help:

- With the feedback controllers, whose integral terms depend on time, the
  solver has to follow the controller dynamics: its timesteps stay at
  10³–10⁴/Ω_ci and it is slower than CVODE, which works for these cases.
- Helium atoms (see above).
- R = 0.99 with any drift, where no steady state was found.
