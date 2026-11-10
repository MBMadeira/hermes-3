#!/usr/bin/env python3
"""Analytic solution of the simple SOL drift model of Kriete et al.

D.M. Kriete et al., Nucl. Fusion 63 (2023) 026022, section 4.2

Steady state, isothermal, constant pitch Theta and binormal ExB velocity v_E,
uniform particle source between the X-points (|y| < y_x) and none beyond:

    d/dy [ n (Theta v_par + v_E) ] = S                         (Kriete eq. 5)
    m_i d/dy [ n (v_par^2 + v_par v_E / Theta) ] = - dp/dy     (Kriete eq. 7)

with p = n m_i c_s^2 and the Bohm-Chodura condition at the targets

    v_par + v_E / Theta = +/- c_s

Normalising with gamma = v_E / (2 c_s Theta), M = v_par / c_s, n* = n / n_max
and y* = y / y_x, the two conservation laws integrate to

    n* (M + 2 gamma)             = F(y*)
    n* (M^2 + 2 gamma M + 1)     = 1 - gamma^2

with F = (gamma + y*)/2 between the X-points and (gamma +/- 1)/2 beyond them.
Writing a = M + 2 gamma (the total velocity normalised to c_s),

    F a^2 - B a + F = 0,   B = 2 gamma F + 1 - gamma^2

whose subsonic root (|a| <= 1) is a = 2F / (B + sqrt(B^2 - 4F^2)), giving

    n* = (B + sqrt(B^2 - 4 F^2)) / 2,   M = a - 2 gamma

The peak density is n_max = 2 S y_x / (Theta c_s) = 2 S s_x / c_s with
s_x = y_x / Theta the parallel length from the mid-point to an X-point.

Derived quantities (valid for |gamma| <= 1/2):
    n_upper / n_lower = (1 + gamma) / (1 - gamma)        (Kriete eq. 10)
    y*_stagnation     = 3 gamma - 4 gamma^3              (Kriete eq. 11)
    y*(v_y = 0)       = -gamma                           (Kriete eq. 12)
"""

import numpy as np


def flux_function(ystar, gamma):
    """Normalised poloidal particle flux F = n* (M + 2 gamma)."""
    ystar = np.asarray(ystar, dtype=float)
    return 0.5 * (gamma + np.clip(ystar, -1.0, 1.0))


def solution(ystar, gamma):
    """Return (nstar, M) at normalised poloidal positions ystar."""
    F = flux_function(ystar, gamma)
    B = 2.0 * gamma * F + 1.0 - gamma**2
    root = np.sqrt(np.maximum(B**2 - 4.0 * F**2, 0.0))
    nstar = 0.5 * (B + root)
    a = 2.0 * F / (B + root)
    return nstar, a - 2.0 * gamma


def target_density_ratio(gamma):
    """n_upper / n_lower (Kriete eq. 10)."""
    return (1.0 + gamma) / (1.0 - gamma)


def stagnation_point(gamma):
    """y* where v_par = 0 (Kriete eq. 11)."""
    return 3.0 * gamma - 4.0 * gamma**3


def zero_poloidal_flux_point(gamma):
    """y* where v_y = Theta v_par + v_E = 0 (Kriete eq. 12)."""
    return -gamma


if __name__ == "__main__":
    # Self-checks of the analytic solution against Kriete's closed forms
    for gamma in [0.0, 0.1, 0.25, -0.25, 0.4, 0.5]:
        n_lo, M_lo = solution(-2.0, gamma)
        n_up, M_up = solution(2.0, gamma)
        assert np.isclose(M_up, 1 - 2 * gamma) and np.isclose(M_lo, -1 - 2 * gamma)
        assert np.isclose(n_up / n_lo, target_density_ratio(gamma))
        ys = stagnation_point(gamma)
        assert abs(solution(ys, gamma)[1]) < 1e-12
        y0 = zero_poloidal_flux_point(gamma)
        n0, M0 = solution(y0, gamma)
        assert abs(M0 + 2 * gamma) < 1e-12
        y = np.linspace(-1, 1, 2001)
        n, M = solution(y, gamma)
        assert np.isclose(n.max(), 1.0, atol=1e-6)
        print(
            f"gamma={gamma:+.2f}: n_up/n_lo={n_up / n_lo:.4f}  "
            f"y*_stag={ys:+.4f}  y*(v_y=0)={y0:+.4f}  "
            f"M_lower={M_lo:+.3f}  M_upper={M_up:+.3f}"
        )
    print("Analytic solution consistent with Kriete eqs. (10)-(12)")
