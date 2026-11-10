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
| v_E = 132 m/s (γ ≈ 0.46) | 2.697 | 2.722 | −0.949, +1.033 | 1.013 |
| v_E = 66 m/s, κ × 30 (T flat to 2 eV) | 1.696 | 1.699 | −0.960, +0.978 | 1.010 |

γ = v_E / (2 c_s Θ) with c_s at the mid-point. Particle balance closes to
better than 10⁻³ in all cases. Because conduction keeps the temperature nearly
uniform here, the target density asymmetry follows Kriete's isothermal result
(1+γ)/(1−γ); it approaches it further as the conductivity is increased.

Other checks: v_E → −v_E mirrors the solution (2.8×10⁻⁵); 1 and 2 processors
agree (1.5×10⁻⁴); ny = 400 gives n_upper/n_lower = 1.629. The case with
binormal anomalous transport (`binormal_stpm` with D = χ = ν = 1 m²/s, listed
before `electron_force_balance`, which uses the electron momentum source) also
reaches a steady state.

In Kriete's geometry with source-free legs (`y_xpoint = 0.4`) the solution is
steady with pressure evolved; the legs become supersonic (u ≈ −1.14, +1.20 c_s)
as the temperature drops towards the targets.

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
| v_E = 132 m/s (γ ≈ 0.47) | 2.705 | 2.749 | −0.935, +1.028 | 1.018 |

Particle balance closes to 3×10⁻⁴; v_E → −v_E mirrors the solution (3×10⁻⁴)
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
