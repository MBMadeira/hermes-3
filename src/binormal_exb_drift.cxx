#include "../include/binormal_exb_drift.hxx"

#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/difops.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/mesh.hxx>
#include <bout/output_bout_types.hxx>

#include <cmath>
#include <string>

BinormalExBSettings readBinormalExBSettings(Options& options, const Options& units) {
  // Reference sound speed [m/s]: rho_s0 * Omega_ci
  const BoutReal Cs0 = get<BoutReal>(units["meters"]) / get<BoutReal>(units["seconds"]);

  const Field3D Theta =
      options["Theta"]
          .doc("Field line pitch inside the island, Theta = sin(theta_i)")
          .withDefault(Field3D(1e-3));

  if (min(Theta, true, "RGN_ALL") <= 0.0) {
    throw BoutException("binormal ExB drift: Theta must be positive everywhere (min {})",
                        min(Theta, true, "RGN_ALL"));
  }

  BinormalExBSettings settings;
  settings.momentum_advection =
      options["momentum_advection"]
          .doc("Advect parallel momentum with the binormal drift?")
          .withDefault<bool>(true);
  settings.energy_advection =
      options["energy_advection"]
          .doc("Advect internal energy (3/2 p) with the binormal drift?")
          .withDefault<bool>(true);

  const std::string mode =
      options["vE_mode"]
          .doc("How the binormal ExB velocity is set: 'constant' (input vE) or "
               "'sheath' (from the sheath potential at the targets)")
          .withDefault<std::string>("constant");

  if (mode == "constant") {
    settings.mode = BinormalExBSettings::Mode::constant;
    const Field3D vE =
        options["vE"]
            .doc("Binormal ExB drift velocity [m/s]. Positive is towards increasing y "
                 "(the upper target)")
            .withDefault(Field3D(0.0));
    settings.W = vE / (Theta * Cs0);
    return settings;
  }
  if (mode != "sheath") {
    throw BoutException("binormal ExB drift: vE_mode must be 'constant' or 'sheath', "
                        "not '{:s}'",
                        mode);
  }

  settings.mode = BinormalExBSettings::Mode::sheath;
  settings.inv_Theta_Cs0 = 1. / (Theta * Cs0);
  settings.W = 0.0 * Theta; // Set in update()
  settings.Tnorm = get<BoutReal>(units["eV"]);

  const BoutReal lambda =
      options["potential_decay_length"]
          .doc("Radial decay length of the potential [m]: E_r = phi_sheath / lambda")
          .as<BoutReal>();
  if (!(lambda > 0.0)) {
    throw BoutException("binormal ExB drift: potential_decay_length must be > 0 ({:e})",
                        lambda);
  }
  const BoutReal direction =
      options["vE_direction"]
          .doc("Direction of v_E: +1 towards increasing y (the upper target), -1 "
               "towards the lower target")
          .withDefault<BoutReal>(1.0);
  if (std::abs(std::abs(direction) - 1.0) > 1e-12) {
    throw BoutException("binormal ExB drift: vE_direction must be +1 or -1 ({:e})",
                        direction);
  }
  const BoutReal B = options["B"]
                         .doc("Magnetic field [T] in v_E = phi / (B lambda)")
                         .withDefault<BoutReal>(get<BoutReal>(units["Tesla"]));
  settings.vE_per_volt = direction / (B * lambda);

  const std::string reference =
      options["sheath_potential"]
          .doc("Sheath potential that sets v_E: 'mean' of the two targets, 'lower' "
               "or 'upper'")
          .withDefault<std::string>("mean");
  if (reference == "mean") {
    settings.reference = BinormalExBSettings::Reference::mean;
  } else if (reference == "lower") {
    settings.reference = BinormalExBSettings::Reference::lower;
  } else if (reference == "upper") {
    settings.reference = BinormalExBSettings::Reference::upper;
  } else {
    throw BoutException("binormal ExB drift: sheath_potential must be 'mean', 'lower' "
                        "or 'upper', not '{:s}'",
                        reference);
  }

  settings.secondary_electron_coef =
      options["secondary_electron_coef"]
          .doc("Secondary electron emission coefficient in the sheath potential")
          .withDefault<BoutReal>(0.0);
  settings.ion_polytropic =
      options["sheath_ion_polytropic"]
          .doc("Ion polytropic coefficient in the Bohm speed for the sheath potential")
          .withDefault<BoutReal>(1.0);
  return settings;
}

std::vector<Permissions::VarRights> BinormalExBSettings::statePermissions() const {
  if (mode != Mode::sheath) {
    return {};
  }
  return {readOnly("species:e:density", Regions::Interior),
          readOnly("species:e:temperature", Regions::Interior),
          readOnly("species:{ions}:density", Regions::Interior),
          readOnly("species:{ions}:temperature", Regions::Interior),
          readOnly("species:{ions}:AA"),
          readOnly("species:{ions}:charge")};
}

void BinormalExBSettings::update(GuardedOptions& state) {
  if (mode != Mode::sheath) {
    return;
  }
  Mesh* mesh = W.getMesh();
  GuardedOptions allspecies = state["species"];
  GuardedOptions electrons = allspecies["e"];
  const Field3D Ne = getNoBoundary<Field3D>(electrons["density"]);
  const Field3D Te = getNoBoundary<Field3D>(electrons["temperature"]);

  // Ion species: density, temperature, mass and charge
  struct Ion {
    Field3D N, T;
    BoutReal AA, Z;
  };
  std::vector<Ion> ions;
  for (auto& kv : allspecies.getChildren()) {
    if (identifySpeciesType(kv.first) != SpeciesType::ion) {
      continue;
    }
    GuardedOptions species = allspecies[kv.first];
    ions.push_back({getNoBoundary<Field3D>(species["density"]),
                    getNoBoundary<Field3D>(species["temperature"]),
                    get<BoutReal>(species["AA"]), get<BoutReal>(species["charge"])});
  }

  const BoutReal Me = SI::Me / SI::Mp;

  // Zero-current sheath potential in the cell (x, y, z), as in
  // sheath_boundary_simple with the Bohm ion flux
  auto sheath_potential = [&](int jx, int jy, int jz) {
    BoutReal ion_sum = 0.0;
    for (const auto& ion : ions) {
      const BoutReal ti = floor(ion.T(jx, jy, jz), 1e-5);
      const BoutReal te = floor(Te(jx, jy, jz), 1e-5);
      ion_sum +=
          ion.Z * ion.N(jx, jy, jz) * sqrt((ion_polytropic * ti + ion.Z * te) / ion.AA);
    }
    const BoutReal te = floor(Te(jx, jy, jz), 1e-5);
    return te
           * log(sqrt(te / (Me * TWOPI)) * (1. - secondary_electron_coef)
                 * floor(Ne(jx, jy, jz), 1e-5) / floor(ion_sum, 1e-5));
  };

  // Sum of potentials and number of cells next to each target
  BoutReal local[4] = {0.0, 0.0, 0.0, 0.0};
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      local[0] += sheath_potential(r.ind, mesh->ystart, jz);
      local[1] += 1.0;
    }
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < mesh->LocalNz; jz++) {
      local[2] += sheath_potential(r.ind, mesh->yend, jz);
      local[3] += 1.0;
    }
  }
  BoutReal global[4];
  MPI_Allreduce(local, global, 4, MPI_DOUBLE, MPI_SUM, BoutComm::get());
  if ((global[1] == 0.0) or (global[3] == 0.0)) {
    throw BoutException("binormal ExB drift: vE_mode = sheath needs targets at both "
                        "y boundaries");
  }
  phi_lower = Tnorm * global[0] / global[1]; // [V]
  phi_upper = Tnorm * global[2] / global[3];

  BoutReal phi = 0.5 * (phi_lower + phi_upper);
  if (reference == Reference::lower) {
    phi = phi_lower;
  } else if (reference == Reference::upper) {
    phi = phi_upper;
  }
  vE = vE_per_volt * phi;

  // Keep W's location and y direction (sheath components store it aligned).
  // Sheath mode is for 1D models, where aligned and non-aligned indices are
  // the same
  BOUT_FOR(i, W.getRegion("RGN_ALL")) { W[i] = vE * inv_Theta_Cs0[i]; }
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

  // With v_E set by the sheath potential, the target plasma is read too. Ion
  // density is already read in the whole domain.
  for (const auto& permission : settings.statePermissions()) {
    if (permission.name != "species:{ions}:density") {
      setPermissions(permission);
    }
  }

  diagnose = options["diagnose"].doc("Output diagnostics?").withDefault<bool>(false);
}

void BinormalExBDrift::transform_impl(GuardedOptions& state) {
  if (settings.selfConsistent()) {
    // v_E from the sheath potentials of the current state
    settings.update(state);
    abs_W = abs(settings.W);
  }

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

  if (settings.selfConsistent()) {
    // Time-dependent v_E set by the sheath potential
    set_with_attrs(state[std::string("w_") + name], settings.W,
                   {{"time_dimension", "t"},
                    {"units", "m / s"},
                    {"conversion", Cs0},
                    {"standard_name", "effective parallel drift velocity"},
                    {"long_name", "Binormal ExB velocity divided by field line pitch"},
                    {"source", "binormal_exb_drift"}});
    set_with_attrs(state[std::string("vE_") + name], settings.vE,
                   {{"time_dimension", "t"},
                    {"units", "m / s"},
                    {"conversion", 1.0},
                    {"standard_name", "binormal ExB velocity"},
                    {"long_name", "Binormal ExB velocity from the sheath potential"},
                    {"source", "binormal_exb_drift"}});
    set_with_attrs(state[std::string("phi_sheath_lower_") + name], settings.phi_lower,
                   {{"time_dimension", "t"},
                    {"units", "V"},
                    {"conversion", 1.0},
                    {"standard_name", "sheath potential"},
                    {"long_name", "Sheath potential at the lower target"},
                    {"source", "binormal_exb_drift"}});
    set_with_attrs(state[std::string("phi_sheath_upper_") + name], settings.phi_upper,
                   {{"time_dimension", "t"},
                    {"units", "V"},
                    {"conversion", 1.0},
                    {"standard_name", "sheath potential"},
                    {"long_name", "Sheath potential at the upper target"},
                    {"source", "binormal_exb_drift"}});
  } else {
    set_with_attrs(state[std::string("w_") + name], settings.W,
                   {{"units", "m / s"},
                    {"conversion", Cs0},
                    {"standard_name", "effective parallel drift velocity"},
                    {"long_name", "Binormal ExB velocity divided by field line pitch"},
                    {"source", "binormal_exb_drift"}});
  }

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
