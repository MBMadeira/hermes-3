# Binormal E×B drift in a 1D island divertor SOL

Benchmark of the `binormal_exb_drift` component, together with the generalised
Bohm condition in `sheath_boundary_simple` (option `exb_drift`), against the
simple SOL drift model of D.M. Kriete et al., Nucl. Fusion 63 (2023) 026022,
section 4.2.

For the same case with evolving ion and electron pressure, see
`../kriete-exb-drift-pressure`.

The model is isothermal, with ion density and parallel momentum evolved, a
constant field-line pitch Θ = 10⁻³, a constant binormal E×B velocity v_E, a
uniform particle source and sheaths at both targets. The drift strength is
γ = v_E / (2 c_s Θ); positive γ is a drift towards the upper target.

## Files

| File | Purpose |
| --- | --- |
| `BOUT.inp` | Input (γ = 0.25, ny = 200, Θ = 10⁻³, L∥ = 1000 m, Tᵢ = Tₑ = 90 eV) |
| `analytic.py` | Kriete's analytic steady state; `python3 analytic.py` runs self-checks against Kriete eqs. (10)–(12) |
| `compare.py` | Compares Hermes-3 runs with the analytic solution and plots profiles |
| `run_benchmark.sh` | Runs the γ scan, resolution scan, island-poloidal grid and Kriete's leg geometry, then compares |

```bash
./run_benchmark.sh /path/to/hermes-3 runs
```

## Results (γ = 0.25)

Analytic values: n_upper/n_lower = (1+γ)/(1−γ) = 1.6667,
stagnation point y*_stag = 3γ − 4γ³ = 0.6875, zero poloidal flux at y* = −γ = −0.25
(y* = y / y_X-point).

| ny | n_upper / n_lower | Stagnation point y* | v_y = 0 at y* | max \|Δn/n_max\| |
| --- | --- | --- | --- | --- |
| 100 | 1.6175 | 0.6635 | −0.2439 | 0.075 |
| 200 | 1.6342 | 0.6745 | −0.2467 | 0.053 |
| 400 | 1.6452 | 0.6803 | −0.2482 | 0.038 |
| 800 | 1.6525 | 0.6835 | −0.2490 | 0.027 |
| analytic | 1.6667 | 0.6875 | −0.2500 | |

Errors decrease steadily with resolution. Convergence is slow because the
analytic solution has a square-root singularity at the sonic points (the
targets); the profile errors are the same with and without the drift (at ny = 200,
max |Δn/n_max| = 0.053 for both γ = 0 and γ = 0.25).

Other checks:

- **Island-poloidal grid.** Running on a grid in y with g²² = Θ² (so that
  `Div_par(F) = d/dy(Θ F)`) gives the same solution as the parallel-length grid
  to within 4×10⁻⁵.
- **MPI.** 1 and 2 processors agree to within 4×10⁻⁵.
- **Bohm condition matters.** With the drift but the standard Bohm condition
  (v∥ = ±c_s, i.e. without `exb_drift` in the sheath), the solution is wrong:
  n_upper/n_lower = 0.71 instead of 1.67, with a subsonic lower target
  (u = −0.49 c_s) and a supersonic upper target (u = 1.39 c_s).
- **Kriete's geometry with source-free legs** (`y_xpoint = 0.4`). The solution
  between the X-points matches the analytic one closely (stagnation point 0.6875,
  v_y = 0 at −0.250 at ny = 400), but the flow in the legs is exactly sonic,
  which is neutrally stable. The solution keeps drifting slowly instead of
  settling: at ny = 400 the legs go from |u| = 1.01 c_s and n_upper/n_lower = 1.69
  after 50 outputs to |u| = 1.05 c_s and 1.73 after 200 outputs. The symmetric
  case without drift does not settle either. For this reason the default input
  places the X-points at the targets.
