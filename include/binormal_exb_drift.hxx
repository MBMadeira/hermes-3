#pragma once
#ifndef BINORMAL_EXB_DRIFT_H
#define BINORMAL_EXB_DRIFT_H

#include <map>
#include <string>

#include <bout/field3d.hxx>
#include <bout/options.hxx>

#include "component.hxx"

/// Effective parallel velocity of a binormal ExB drift in a 1D island SOL model
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
///
/// Reads, from the given options section:
///   - vE     Binormal ExB velocity [m/s]. Positive is towards increasing y
///            (towards the upper target). May be a function of y.
///   - Theta  Field line pitch inside the island. May be a function of y.
///
/// Returns w normalised to the reference sound speed Cs0.
///
/// This is a free function so that sheath boundary components can compute
/// exactly the same velocity as the drift component (they cannot read it from
/// the state without creating a circular dependency).
Field3D binormalExBParallelVelocity(Options& options, const Options& units);

/// Binormal (island-poloidal) ExB drift for 1D stellarator SOL models
///
/// Adds to the ion density and parallel momentum equations
///
///     dn/dt     += - Div_par( n w )
///     d(nv)/dt  += - Div_par( m n v_par w )
///
/// with w = v_E / Theta (see binormalExBParallelVelocity). The fluxes are
/// upwinded (Lax flux with wave speed |w|, MC-limited reconstruction), and at
/// sheath boundaries use the boundary (guard-cell average) values set by the
/// sheath component. The total particle flux through a target is then
/// n (v_par + w), which is what the generalised Bohm condition in
/// sheath_boundary_simple (option exb_drift) imposes to be +/- n c_s.
///
/// The operators are written in terms of Div_par, so the component works
/// both on a grid in parallel length s_par (g_22 = 1, constant Theta) and on
/// a grid in island-poloidal length s_theta with g_22 = 1/Theta^2, where
/// Div_par(n w) = d/ds_theta (n v_E) for any Theta(s_theta).
///
/// Energy terms are not yet included: use with isothermal species, or note
/// that the ExB advection of internal energy is missing.
struct BinormalExBDrift : public NamedComponent<BinormalExBDrift> {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - vE                  Binormal ExB velocity [m/s]
  ///   - Theta               Field line pitch
  ///   - momentum_advection  Include the momentum advection term? Default true
  ///   - diagnose            Save diagnostics? Default false
  ///
  BinormalExBDrift(std::string name, Options& alloptions, Solver* solver);

  void outputVars(Options& state) override;

  static constexpr auto type = "binormal_exb_drift";

private:
  Field3D W;               ///< Effective parallel drift velocity v_E / Theta, normalised
  Field3D abs_W;           ///< |W|, wave speed used in the Lax flux
  bool momentum_advection; ///< Include the momentum advection term?
  bool diagnose;           ///< Save diagnostics?

  /// Diagnostics: particle source due to the drift, and particle flow
  /// through lower y cell faces, for each species
  std::map<std::string, Field3D> particle_source, particle_flow_ylow;

  /// Sets
  /// - species
  ///   - <ion species>
  ///     - density_source
  ///     - momentum_source
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<BinormalExBDrift> registercomponentbinormalexbdrift;
}

#endif // BINORMAL_EXB_DRIFT_H
