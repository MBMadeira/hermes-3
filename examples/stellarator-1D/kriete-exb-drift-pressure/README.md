# Binormal E×B drift with evolving pressure

Extends the isothermal benchmark in `../kriete-exb-drift` (Kriete et al., Nucl.
Fusion 63 (2023) 026022) by evolving the ion and electron pressure, with
Braginskii parallel conduction, collisions, friction and heat exchange.
Uniform particle and power sources, sheaths at both targets, Θ = 10⁻³,
L∥ = 1000 m.

The drift enters the pressure equation of each species as

    d(3/2 p)/dt = ... − Div_par(5/2 p w) + w Grad_par(p),     w = v_E / Θ

which equals −Div_par(3/2 p w) in the interior (the E×B drift is
incompressible). It is written in the same form as the parallel advection in
`evolve_pressure`, so that the energy flux through a target is the enthalpy
5/2 p w carried by the particles. In `sheath_boundary_simple` (option
`exb_drift`) the Bohm condition, the sheath potential and the sheath heat
fluxes use the total velocity u = v∥ + w, and the ion heat flux uses the actual
parallel kinetic energy at the sheath, ½ m v∥² with v∥ = u − w.

## Running

```bash
hermes-3 -d .
python3 analyse.py . --plot profiles.png
```

`analyse.py` prints the target temperatures and densities, u/c_s at both
targets, and the particle and energy balance. Energy enters through the
power source and the E×B work ∫ w ∂p ds = w (p_upper − p_lower), and leaves
through the sheaths and as the kinetic energy given to source particles
(injected at rest; Hermes-3 does not return this to heat).

## Results (ny = 200, 50 outputs)

| Case | n_upper / n_lower | Kriete (isothermal) | u/c_s at targets | Energy balance |
| --- | --- | --- | --- | --- |
| v_E = 0 | 1.000 | 1.000 | −0.978, +0.978 | 1.006 |
| v_E = 66 m/s (γ ≈ 0.23) | 1.616 | 1.612 | −0.965, +0.989 | 1.008 |
| v_E = −66 m/s | 0.619 | 0.620 | −0.989, +0.965 | 1.008 |
| v_E = 132 m/s (γ ≈ 0.46) | 2.554 | 2.724 | −0.948, +1.144 | 1.013 |
| v_E = 66 m/s, κ × 30 (T flat to 2 eV) | 1.696 | 1.699 | −0.960, +0.978 | 1.010 |

γ = v_E / (2 c_s Θ) with c_s at the mid-point. Particle balance closes to
better than 10⁻³ in all cases. Because conduction keeps the temperature nearly
uniform here, the target density asymmetry follows Kriete's isothermal result
(1+γ)/(1−γ); it approaches it further as the conductivity is increased. The
exception is v_E = 132 m/s: there w = v_E/Θ = 132 km/s exceeds c_s ≈ 127 km/s at
the upper target, so the drift alone carries the plasma into that target
supersonically (u = 1.14 c_s), and the ratio is 6% below Kriete's.

Other checks: v_E → −v_E mirrors the solution (2.8×10⁻⁵); 1 and 2 processors
agree (1.5×10⁻⁴); ny = 400 gives n_upper/n_lower = 1.629. The case with
binormal anomalous transport (`binormal_stpm` with D = χ = ν = 1 m²/s, listed
before `electron_force_balance`, which uses the electron momentum source) also
reaches a steady state.

In Kriete's geometry with source-free legs (`y_xpoint = 0.4`) the solution is
steady with pressure evolved; the legs become supersonic (u ≈ −1.43, +1.58 c_s)
as the temperature drops towards the targets. When u > c_s the ion heat flux
uses the actual parallel kinetic energy ½ m v∥² in place of the ½ m c_s²
assumed in γ_i (an earlier version subtracted ½ m u², which is only correct
at u = c_s, and gave u ≈ −1.14, +1.20 c_s here).

The energy balance error (≲ 1.5%) is of the same size as without the drift
and comes mainly from approximating the face pressures in the E×B work by the
values in the cells next to the targets.

## With the full sheath model (`sheath_boundary`)

The `exb_drift` option is also available in `sheath_boundary` (Tskhakaya
ion speed and heat transmission, sheath potential, secondary emission). To use
it, replace `sheath_boundary_simple` by `sheath_boundary` in
`hermes:components` and rename the `[sheath_boundary_simple]` section; both
options (`exb_drift`, `diagnose`) are the same. With
`analyse.py --bohm-adiabatic 1.6667` the Bohm check uses the ion speed of
that component, C_i² = (5/3 T_i + T_e)/m_i.

| Case (ny = 200) | n_upper / n_lower | Kriete (isothermal) | u/C_i at targets | Energy balance |
| --- | --- | --- | --- | --- |
| v_E = 0 | 1.000 | 1.000 | −0.969, +0.969 | 1.012 |
| v_E = 66 m/s (γ ≈ 0.24) | 1.607 | 1.617 | −0.954, +0.982 | 1.014 |
| v_E = −66 m/s | 0.623 | 0.618 | −0.982, +0.954 | 1.014 |
| v_E = 132 m/s (γ ≈ 0.46) | 2.747 | 2.733 | −0.935, +1.018 | 1.018 |

Particle balance closes to 3×10⁻⁴ (1% at v_E = 132 m/s, which is not fully
steady at 50 outputs); v_E → −v_E mirrors the solution (3×10⁻⁴)
and 1 and 2 processors agree (7×10⁻⁵). The energy balance is ~1% higher than
with the simple sheath already without drift: the total sheath loss diagnosed by
`sheath_boundary` approximates the advected part with sheath-entrance values.
Without drift the boundary conditions are unchanged (identical right-hand side).

At ny = 400 the runs with this sheath do not settle completely: the density
near the targets fluctuates by ~10⁻³ between outputs without drift, and by up
to ~10⁻² with v_E = 66 m/s (near the target the drift points away from). The
simple sheath at the same resolution is steady. In the isothermal benchmark
(`../kriete-exb-drift`) this sheath imposes C_i = 1.15 c_s at the targets
(5/3 T_i in the Bohm speed), which gives a supersonic layer of two cells there;
the solution between the targets is close to Kriete's (n_upper/n_lower = 1.632,
stagnation point 0.668 at ny = 200).

## v_E from the sheath potential

With `binormal_exb_drift:vE_mode = sheath`, v_E = ±φ/(B λ_φ) is set by the
radial electric field φ/λ_φ of the sheath potential: φ is the zero-current
sheath potential (as in `sheath_boundary_simple`, with the Bohm ion flux) in
the cells next to the targets, averaged over the two targets, and λ_φ is the
radial decay length of the potential (`potential_decay_length`). v_E is uniform
along the field line and follows the target plasma in time. See the commented
lines in `[binormal_exb_drift]` in `BOUT.inp`.

Here λ_φ = 1.136 m was chosen so that the steady state with constant
v_E = 66 m/s gives v_E = 66 m/s (φ ≈ 187 V at T_e ≈ 76 eV, B = 2.5 T). With
that λ_φ, the sheath-mode run reaches the same steady state (profiles within
5×10⁻⁵, v_E = 66.00 m/s after a few outputs). Changing the power, at fixed
λ_φ (ny = 200, simple sheath):

| Power (both species) | v_E mode | v_E [m/s] | γ | T_e lower / upper [eV] | n_upper / n_lower |
| --- | --- | --- | --- | --- | --- |
| 0.6×10⁸ W/m² | constant | 66 | 0.32 | 35.6 / 33.4 | 1.37 |
| 0.6×10⁸ W/m² | sheath | 29.7 | 0.14 | 35.6 / 33.5 | 1.05 |
| 1.2×10⁸ W/m² | constant | 66 | 0.23 | 77.0 / 75.4 | 1.616 |
| 1.2×10⁸ W/m² | sheath | 66.0 | 0.23 | 77.0 / 75.4 | 1.616 |
| 2.4×10⁸ W/m² | constant | 66 | 0.18 | 152.6 / 147.8 | 1.44 |
| 2.4×10⁸ W/m² | sheath | 133.2 | 0.35 | 159.1 / 149.5 | 2.13 |

With constant v_E, γ = v_E/(2 c_s Θ) grows as the plasma cools; with v_E from
the sheath potential, v_E ∝ T_e at the targets and γ ∝ √T_e, so the drift
asymmetry grows with power instead. At half power the upper target is
supersonic (u/c_s ≈ 1.3–1.5 in the last cell) and the asymmetry is below
Kriete's isothermal estimate in both modes.

The decay length needed for v_E ≈ 66 m/s with Θ = 10⁻³ is about 1 m; for a
decay length of a few cm, v_E from this formula is tens of times larger and
γ ≫ 1. The values of λ_φ and Θ (and B) therefore need to be chosen together.
