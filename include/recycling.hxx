#pragma once
#ifndef RECYCLING_H
#define RECYCLING_H

#include "binormal_exb_drift.hxx"
#include "component.hxx"

/// Convert fluxes of species at boundaries
///
/// Since this must be calculated after boundary fluxes (e.g. sheath),
/// it is included as a top-level component
///
struct Recycling : public NamedComponent<Recycling> {

  /// Inputs
  ///
  ///   - <name>
  ///     - species    A comma-separated list of species to recycle
  ///   - <species>
  ///     - recycle_as  The species to recycle into
  ///     - recycle_multiplier   The recycled flux multiplier, between 0 and 1
  ///     - recycle_energy       The energy of the recycled particles [eV]
  ///
  ///   - <name>
  ///     - exb_drift  Name of a binormal_exb_drift component section. If set, the
  ///                  target particle flux that is recycled is n (v_par + v_E/Theta),
  ///                  the total flux into the target including the binormal ExB drift,
  ///                  and fast recycling uses gamma_i T n (v_par + v_E/Theta).
  ///                  The sheath should use the same drift (sheath_boundary_simple:
  ///                  exb_drift or sheath_boundary:exb_drift).
  ///
  Recycling(std::string name, Options& alloptions, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "recycling";

private:
  struct RecycleChannel {
    std::string from; ///< The species name to recycle
    std::string to;   ///< Species to recycle to

    /// Flux multiplier (recycling fraction).
    /// Combination of recycling fraction and species change e.g h+ -> h2 results in 0.5
    /// multiplier
    BoutReal target_multiplier, sol_multiplier, pfr_multiplier, pump_multiplier;
    BoutReal target_energy, sol_energy,
        pfr_energy; ///< Energy of recycled particle (normalised to Tnorm)
    BoutReal target_fast_recycle_fraction, pfr_fast_recycle_fraction,
        sol_fast_recycle_fraction; ///< Fraction of ions undergoing fast reflection
    BoutReal target_fast_recycle_energy_factor, sol_fast_recycle_energy_factor,
        pfr_fast_recycle_energy_factor; ///< Fraction of energy retained by fast recycled
                                        ///< neutrals

    // Recycling particle and energy sources for the different sources of recycling
    // These sources are per-channel and added to the `to` species
    /// Recycling density source for target recycling
    Field3D target_recycle_density_source = 0.0;
    /// Recycling energy source for target recycling
    Field3D target_recycle_energy_source = 0.0;
    /// Recycling density source for pfr + sol recycling
    Field3D wall_recycle_density_source = 0.0;
    /// Recycling energy source for pfr + sol recycling
    Field3D wall_recycle_energy_source = 0.0;
    /// Recycling density source for pump recycling
    Field3D pump_density_source = 0.0;
    /// Recycling energy source for pump recycling
    Field3D pump_energy_source = 0.0;
  };

  std::vector<RecycleChannel> channels; // Recycling channels

  bool target_recycle{false}, sol_recycle{false}, pfr_recycle{false},
      neutral_pump{false}; ///< Flags for enabling recycling in different regions
  bool diagnose;           ///< Save additional post-processing variables?
  bool has_sheath_boundary_simple{
      false}; ///< Whether sheath_boundary_simple component is available

  BoutReal density_floor,
      pressure_floor; ///< minimum values for Nn, Pn to avoid divide by zero

  BoutReal gamma_i; /// Sheath heat transmission coefficient

  /// Binormal ExB drift (see binormal_exb_drift.hxx). The drift's effective parallel
  /// velocity W = v_E/Theta is added to the ion parallel velocity in the target
  /// flux. Zero if no drift.
  bool use_exb_drift{false};
  BinormalExBSettings exb;

  Field3D density_source,
      energy_source;          ///< Recycling particle and energy sources for all locations
  Field3D energy_flow_xlow;   ///< Cell edge fluxes used for calculating
                              ///<
                              ///< fast recycling energy source
  Field3D particle_flow_xlow; ///< Radial wall particle fluxes for recycling calc. No need
                              ///<
                              ///< to get poloidal from here, it's calculated from sheath
                              ///<
                              ///< velocity

  Field2D is_pump; ///< 1 = pump, 0 = no pump. Provided by user in
                   ///<
                   ///< grid file.
  /// Inputs
  ///
  /// - species
  ///   - <from species>
  ///     - density
  ///     - velocity
  ///     - temperature
  ///   - <to species>
  ///     - AA
  ///     - density
  ///     - pressure
  ///     - temperature
  ///
  /// Outputs
  ///
  /// - species
  ///   - <species>
  ///     - density_source
  ///     - energy_source
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<Recycling> registercomponentrecycling;
}

#endif // RECYCLING_H
