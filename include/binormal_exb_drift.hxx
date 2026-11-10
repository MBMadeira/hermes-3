#pragma once
#ifndef BINORMAL_EXB_DRIFT_H
#define BINORMAL_EXB_DRIFT_H

#include <map>
#include <string>

#include <bout/field3d.hxx>
#include <bout/options.hxx>

#include "component.hxx"

/// Settings of a binormal ExB drift in a 1D island SOL model
///
/// In a stellarator island the field line pitch inside the island is
/// Theta = sin(theta_i) ~ 1e-3. With helical symmetry (d/d zeta = 0)
/// a binormal flux F_perp contributes to the 1D divergence like a
/// parallel flux F_perp / Theta. The binormal ExB velocity v_E therefore
/// acts as an additional parallel advection velocity
///
///     w = v_E / Theta
///
/// See Kriete et al., Nucl. Fusion 63 (2023) 026022, eqs (5) and (7).
struct BinormalExBSettings {
  Field3D W;                      ///< w = v_E / Theta, normalised to Cs0
  bool momentum_advection{false}; ///< Advect parallel momentum with the drift?
  bool energy_advection{false};   ///< Advect internal energy with the drift?

  /// Energy flux through a target carried by the binormal_exb_drift component,
  /// for a species with temperature T, density n, parallel velocity vpar and
  /// mass `mass` at the sheath entrance, where the drift velocity is w:
  /// enthalpy 5/2 n T w (if energy is advected) plus kinetic energy
  /// 1/2 m n v_par^2 w (if momentum is advected).
  ///
  /// Sheath components subtract this from their heat flux, so that the total
  /// energy leaving through the target is the sheath heat flux.
  BoutReal targetEnergyFlux(BoutReal T, BoutReal n, BoutReal vpar, BoutReal mass,
                            BoutReal w) const {
    BoutReal flux = 0.0;
    if (energy_advection) {
      flux += 2.5 * T * n * w;
    }
    if (momentum_advection) {
      flux += 0.5 * mass * vpar * vpar * n * w;
    }
    return flux;
  }
};

/// Read the binormal ExB drift settings from an options section
///
/// Reads:
///   - vE                  Binormal ExB velocity [m/s]. Positive is towards
///                         increasing y (the upper target). May depend on y.
///   - Theta               Field line pitch inside the island. May depend on y.
///   - momentum_advection  Advect parallel momentum? Default true
///   - energy_advection    Advect internal energy? Default true
///
/// This is a free function so that sheath boundary components use exactly the
/// same settings as the drift component (they cannot read them from the state
/// without creating a circular dependency).
BinormalExBSettings readBinormalExBSettings(Options& options, const Options& units);

/// Binormal (island-poloidal) ExB drift for 1D stellarator SOL models
///
/// Adds, with w = v_E / Theta,
///
///     dn/dt          += - Div_par( n w )                         (ions)
///     d(m n v)/dt    += - Div_par( m n v_par w )                 (ions)
///     d(3/2 p)/dt    += - Div_par( 5/2 p w ) + w Grad_par(p)     (ions and electrons)
///
/// The ExB drift is the same for all charged species and is incompressible,
/// so in the interior the energy terms equal - Div_par(3/2 p w): the work term
/// w Grad(p) cancels part of the enthalpy flux. They are written in the same
/// form as the parallel advection in evolve_pressure so that the energy flux
/// through a target is the enthalpy flux 5/2 p w carried by the particles.
///
/// The fluxes are upwinded (Lax flux with wave speed |w|, MC-limited
/// reconstruction), and at sheath boundaries use the boundary (guard-cell
/// average) values set by the sheath component. The total particle flux
/// through a target is then n (v_par + w), which the generalised Bohm
/// condition in sheath_boundary_simple (option exb_drift) sets to +/- n c_s.
///
/// The operators are written in terms of Div_par, so the component works
/// both on a grid in parallel length s_par (g_22 = 1, constant Theta) and on
/// a grid in island-poloidal length s_theta with g^22 = Theta^2, where
/// Div_par(n w) = d/ds_theta (n v_E) for any Theta(s_theta).
struct BinormalExBDrift : public NamedComponent<BinormalExBDrift> {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - vE                  Binormal ExB velocity [m/s]
  ///   - Theta               Field line pitch
  ///   - momentum_advection  Include the momentum advection term? Default true
  ///   - energy_advection    Include the energy advection term? Default true
  ///   - diagnose            Save diagnostics? Default false
  ///
  BinormalExBDrift(std::string name, Options& alloptions, Solver* solver);

  void outputVars(Options& state) override;

  static constexpr auto type = "binormal_exb_drift";

private:
  BinormalExBSettings settings;
  Field3D abs_W; ///< |w|, wave speed used in the Lax flux
  bool diagnose; ///< Save diagnostics?

  /// Diagnostics for each species: particle and energy sources due to the
  /// drift, and particle and energy flows through lower y cell faces
  std::map<std::string, Field3D> particle_source, particle_flow_ylow, energy_source,
      energy_flow_ylow;

  /// Sets
  /// - species
  ///   - <ion species>
  ///     - density_source
  ///     - momentum_source
  ///     - energy_source
  ///   - e
  ///     - energy_source
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<BinormalExBDrift> registercomponentbinormalexbdrift;
}

#endif // BINORMAL_EXB_DRIFT_H
