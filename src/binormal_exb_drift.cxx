#include "../include/binormal_exb_drift.hxx"

#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/difops.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include <string>

BinormalExBSettings readBinormalExBSettings(Options& options, const Options& units) {
  // Reference sound speed [m/s]: rho_s0 * Omega_ci
  const BoutReal Cs0 = get<BoutReal>(units["meters"]) / get<BoutReal>(units["seconds"]);

  const Field3D vE =
      options["vE"]
          .doc("Binormal ExB drift velocity [m/s]. Positive is towards increasing y "
               "(the upper target)")
          .withDefault(Field3D(0.0));

  const Field3D Theta =
      options["Theta"]
          .doc("Field line pitch inside the island, Theta = sin(theta_i)")
          .withDefault(Field3D(1e-3));

  if (min(Theta, true, "RGN_ALL") <= 0.0) {
    throw BoutException("binormal ExB drift: Theta must be positive everywhere (min {})",
                        min(Theta, true, "RGN_ALL"));
  }

  return {vE / (Theta * Cs0),
          options["momentum_advection"]
              .doc("Advect parallel momentum with the binormal drift?")
              .withDefault<bool>(true),
          options["energy_advection"]
              .doc("Advect internal energy (3/2 p) with the binormal drift?")
              .withDefault<bool>(true)};
}

BinormalExBDrift::BinormalExBDrift(std::string name, Options& alloptions,
                                   [[maybe_unused]] Solver* solver)
    : NamedComponent(name, {
                               // Density, momentum and pressure are read in the whole
                               // domain: the boundary (guard cell) values set by the
                               // sheath determine the drift fluxes through the targets.
                               readIfSet("species:{ions}:density"),
                               readIfSet("species:{ions}:momentum"),
                               readIfSet("species:{charged}:pressure"),
                               readWrite("species:{ions}:density_source"),
                               readWrite("species:{ions}:momentum_source"),
                               readWrite("species:{charged}:energy_source"),
                           }) {
  auto& options = alloptions[name];

  settings = readBinormalExBSettings(options, alloptions["units"]);
  abs_W = abs(settings.W);

  diagnose = options["diagnose"].doc("Output diagnostics?").withDefault<bool>(false);
}

void BinormalExBDrift::transform_impl(GuardedOptions& state) {
  GuardedOptions allspecies = state["species"];
  const Field3D& W = settings.W;

  for (auto& kv : allspecies.getChildren()) {
    const std::string& species_name = kv.first;

    // The drift is the same for all charged species. Ion density and momentum
    // are advected; electron density follows from quasineutrality and electron
    // parallel momentum from force balance. The internal energy of both ions
    // and electrons is advected.
    const SpeciesType species_type = identifySpeciesType(species_name);
    if (species_type == SpeciesType::neutral) {
      continue;
    }
    const bool is_ion = species_type == SpeciesType::ion;

    GuardedOptions species = allspecies[species_name];

    if (is_ion and species.isSet("density")) {
      const Field3D N = GET_VALUE(Field3D, species["density"]);

      Field3D flow_ylow;
      const Field3D div_flux =
          FV::Div_par_mod<hermes::Limiter>(N, W, abs_W, flow_ylow, true);
      subtract(species["density_source"], div_flux);

      if (diagnose) {
        particle_source[species_name] = -div_flux;
        particle_flow_ylow[species_name] = flow_ylow;
      }
    }

    if (is_ion and settings.momentum_advection and species.isSet("momentum")) {
      const Field3D NV = GET_VALUE(Field3D, species["momentum"]);

      Field3D flow_ylow;
      subtract(species["momentum_source"],
               FV::Div_par_mod<hermes::Limiter>(NV, W, abs_W, flow_ylow, true));
    }

    if (settings.energy_advection and species.isSet("pressure")) {
      const Field3D P = GET_VALUE(Field3D, species["pressure"]);

      // Internal energy: d(3/2 p)/dt = ... - Div(5/2 p w) + w Grad(p)
      //
      // This has the same structure as the parallel advection in
      // evolve_pressure: the enthalpy flux 5/2 p w is advected
      // conservatively and the work term w Grad(p) is a volumetric
      // exchange. For incompressible w this equals -Div(3/2 p w) in the
      // interior, but the flux through a target is the enthalpy flux
      // 5/2 p w carried by the particles, consistent with the sheath.
      Field3D flow_ylow;
      const Field3D div_flux =
          2.5 * FV::Div_par_mod<hermes::Limiter>(P, W, abs_W, flow_ylow, true)
          - W * Grad_par(P.asField3DParallel());
      subtract(species["energy_source"], div_flux);

      if (diagnose) {
        energy_source[species_name] = -div_flux;
        energy_flow_ylow[species_name] = 2.5 * flow_ylow;
      }
    }
  }
}

void BinormalExBDrift::outputVars(Options& state) {
  if (not diagnose) {
    return;
  }

  const auto Nnorm = get<BoutReal>(state["Nnorm"]);
  const auto Tnorm = get<BoutReal>(state["Tnorm"]);
  const auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  const auto Cs0 = get<BoutReal>(state["Cs0"]);
  const auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  const std::string& name = objectName();

  set_with_attrs(state[std::string("w_") + name], settings.W,
                 {{"units", "m / s"},
                  {"conversion", Cs0},
                  {"standard_name", "effective parallel drift velocity"},
                  {"long_name", "Binormal ExB velocity divided by field line pitch"},
                  {"source", "binormal_exb_drift"}});

  for (const auto& [species_name, source] : particle_source) {
    set_with_attrs(state[std::string("S") + species_name + std::string("_exb")], source,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "particle source"},
                    {"long_name", species_name + " particle source due to binormal ExB"},
                    {"species", species_name},
                    {"source", "binormal_exb_drift"}});
  }
  for (const auto& [species_name, flow] : particle_flow_ylow) {
    set_with_attrs(state[std::string("F") + species_name + std::string("_exb_ylow")],
                   flow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * rho_s0 * rho_s0 * Nnorm * Omega_ci},
                    {"standard_name", "particle flow"},
                    {"long_name", species_name + " binormal ExB particle flow"},
                    {"species", species_name},
                    {"source", "binormal_exb_drift"}});
  }
  for (const auto& [species_name, source] : energy_source) {
    set_with_attrs(state[std::string("E") + species_name + std::string("_exb")], source,
                   {{"time_dimension", "t"},
                    {"units", "W / m^3"},
                    {"conversion", SI::qe * Tnorm * Nnorm * Omega_ci},
                    {"standard_name", "energy source"},
                    {"long_name", species_name + " energy source due to binormal ExB"},
                    {"species", species_name},
                    {"source", "binormal_exb_drift"}});
  }
  for (const auto& [species_name, flow] : energy_flow_ylow) {
    set_with_attrs(
        state[std::string("EF") + species_name + std::string("_exb_ylow")], flow,
        {{"time_dimension", "t"},
         {"units", "W"},
         {"conversion", rho_s0 * rho_s0 * rho_s0 * SI::qe * Tnorm * Nnorm * Omega_ci},
         {"standard_name", "energy flow"},
         {"long_name", species_name + " binormal ExB energy flow"},
         {"species", species_name},
         {"source", "binormal_exb_drift"}});
  }
}
