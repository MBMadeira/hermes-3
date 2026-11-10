#include "../include/binormal_exb_drift.hxx"

#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <bout/boutexception.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/output_bout_types.hxx>

#include <string>

Field3D binormalExBParallelVelocity(Options& options, const Options& units) {
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

  return vE / (Theta * Cs0);
}

BinormalExBDrift::BinormalExBDrift(std::string name, Options& alloptions,
                                   [[maybe_unused]] Solver* solver)
    : NamedComponent(name, {
                               // Density and momentum are read in the whole domain:
                               // the boundary (guard cell) values set by the sheath
                               // determine the drift flux through the targets.
                               readIfSet("species:{ions}:density"),
                               readIfSet("species:{ions}:momentum"),
                               readWrite("species:{ions}:density_source"),
                               readWrite("species:{ions}:momentum_source"),
                           }) {
  auto& options = alloptions[name];

  W = binormalExBParallelVelocity(options, alloptions["units"]);
  abs_W = abs(W);

  momentum_advection = options["momentum_advection"]
                           .doc("Advect parallel momentum with the binormal drift?")
                           .withDefault<bool>(true);

  diagnose = options["diagnose"].doc("Output diagnostics?").withDefault<bool>(false);
}

void BinormalExBDrift::transform_impl(GuardedOptions& state) {
  GuardedOptions allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    const std::string& species_name = kv.first;

    // The drift is the same for all charged species. Only ion density and
    // momentum are evolved here: electron density follows from quasineutrality
    // and electron parallel momentum from force balance.
    if (identifySpeciesType(species_name) != SpeciesType::ion) {
      continue;
    }

    GuardedOptions species = allspecies[species_name];

    if (species.isSet("density")) {
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

    if (momentum_advection and species.isSet("momentum")) {
      const Field3D NV = GET_VALUE(Field3D, species["momentum"]);

      Field3D flow_ylow;
      subtract(species["momentum_source"],
               FV::Div_par_mod<hermes::Limiter>(NV, W, abs_W, flow_ylow, true));
    }
  }
}

void BinormalExBDrift::outputVars(Options& state) {
  if (not diagnose) {
    return;
  }

  const auto Nnorm = get<BoutReal>(state["Nnorm"]);
  const auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  const auto Cs0 = get<BoutReal>(state["Cs0"]);
  const auto rho_s0 = get<BoutReal>(state["rho_s0"]);
  const std::string& name = objectName();

  set_with_attrs(state[std::string("w_") + name], W,
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
}
