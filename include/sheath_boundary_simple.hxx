#pragma once
#ifndef SHEATH_BOUNDARY_SIMPLE_H
#define SHEATH_BOUNDARY_SIMPLE_H

#include "binormal_exb_drift.hxx"
#include "component.hxx"

/// Boundary condition at the wall in Y
///
/// This is a collective component, because it couples all charged species
///
/// This implements a simple boundary condition, where each species
/// goes to their own sound velocity at the sheath entrance.
///
/// Notes:
///   - It is recommended to use SheathBoundary rather than SheathBoundarySimple;
///     this is here for comparison to that more complete model.
///
struct SheathBoundarySimple : public NamedComponent<SheathBoundarySimple> {
  /// # Input options
  /// - <name>  e.g. "sheath_boundary_simple"
  ///   - lower_y                  Boundary on lower y?
  ///   - upper_y                  Boundary on upper y?
  ///   - gamma_e                  Electron sheath heat transmission coefficient
  ///   - gamma_i                  Ion sheath heat transmission coefficient
  ///   - sheath_ion_polytropic    Ion polytropic coefficient in Bohm sound speed. Default 1.
  ///   - wall_potential           Voltage of the wall [Volts]
  ///   - secondary_electron_coef  Effective secondary electron emission coefficient
  ///   - sin_alpha                Sine of the angle between magnetic field line and wall surface (0 to 1)
  ///   - always_set_phi           Always set phi field? Default is to only modify if already set
  ///   - exb_drift                Name of a binormal_exb_drift component section. If set,
  ///                              the Bohm condition is imposed on the total effective
  ///                              parallel velocity u = v_par + v_E/Theta (see below).
  SheathBoundarySimple(std::string name, Options& options, Solver*);

  void outputVars(Options& state) override;

  static constexpr auto type = "sheath_boundary_simple";

private:
  BoutReal Ge;        // Secondary electron emission coefficient
  BoutReal sin_alpha; // sin of angle between magnetic field and wall.

  BoutReal gamma_e;               ///< Electron sheath heat transmission
  BoutReal gamma_i;               ///< Ion sheath heat transmission
  BoutReal sheath_ion_polytropic; ///< Polytropic coefficient in sheat velocity

  bool lower_y; // Boundary on lower y?
  bool upper_y; // Boundary on upper y?

  bool always_set_phi; ///< Set phi field?

  Field3D wall_potential; ///< Voltage of the wall. Normalised units.

  Field3D hflux_e; // Electron heat flux through sheath
  Field3D phi;     // Phi at sheath
  Field3D ion_sum; // Sum of ion current at sheath

  bool diagnose;       // Save diagnostic variables?
  Options diagnostics; // Options object to store diagnostic fields like a dict

  bool no_flow; ///< No flow speed, only remove energy

  /// Binormal ExB drift settings (see binormal_exb_drift.hxx). The effective
  /// parallel velocity w = v_E / Theta is stored field aligned, and is zero if
  /// no drift is configured.
  bool use_exb_drift{false};
  BinormalExBSettings exb;

  BoutReal density_boundary_mode, pressure_boundary_mode,
      temperature_boundary_mode; ///< BC mode: 0=LimitFree, 1=ExponentialFree, 2=LinearFree

  ///
  /// # Inputs
  /// - species
  ///   - e
  ///     - density
  ///     - temperature
  ///     - pressure    Optional
  ///     - velocity    Optional
  ///     - mass        Optional
  ///   - <ions>  if charge is set (i.e. not neutrals)
  ///     - charge
  ///     - mass
  ///     - density
  ///     - temperature
  ///     - pressure     Optional
  ///     - velocity     Optional. Default 0
  ///     - momentum     Optional. Default mass * density * velocity
  /// - fields
  ///   - phi    Optional. If not set, calculated at boundary (see note below)
  ///
  /// # Outputs
  /// - species
  ///   - e
  ///     - density      Sets boundary
  ///     - temperature  Sets boundary
  ///     - velocity     Sets boundary
  ///     - energy_source
  ///   - <ions>
  ///     - density      Sets boundary
  ///     - temperature  Sets boundary
  ///     - velocity     Sets boundary
  ///     - momentum     Sets boundary
  ///     - energy_source
  /// - fields
  ///   - phi   Sets boundary
  ///
  /// With a binormal ExB drift (option exb_drift), the Bohm-Chodura condition
  /// applies to the velocity normal to the target, which in the 1D island model
  /// is u = v_par + w with w = v_E / Theta (Kriete et al. NF 2023, Stangeby &
  /// Chankin). The ion parallel velocity at the sheath entrance is then
  /// v_par = u - w with |u| >= c_s, so that the total particle flux through the
  /// target n (v_par + w) is the Bohm flux. The drift flux itself is carried by
  /// the binormal_exb_drift component. The sheath potential and the heat fluxes
  /// use the total (normal) flux n u. The total ion energy loss is
  /// gamma_i T n u + 1/2 m (v_par^2 - u^2) n u: gamma_i assumes ions arrive with
  /// kinetic energy 1/2 m u^2, whereas with the drift their parallel velocity is
  /// v_par = u - w (identical without drift). The energy carried through the
  /// target by the drift component (5/2 p w and 1/2 m n v_par^2 w) is subtracted
  /// from the sheath sink. The electron parallel velocity at the sheath is set to
  /// u - w so that J_par = 0.
  ///
  /// If the field phi is set, then this is used in the boundary condition.
  /// If not set, phi at the boundary is calculated and stored in the state.
  /// Note that phi in the domain will not be set, so will be invalid data.
  ///
  ///
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<SheathBoundarySimple> registercomponentsheathboundarysimple;
}

#endif // SHEATH_BOUNDARY_SIMPLE_H
