#pragma once
#ifndef FEEDBACK_MEASUREMENT_H
#define FEEDBACK_MEASUREMENT_H

#include <string>

#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/field3d.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/region.hxx>

/// Where a feedback controller measures the quantity it controls.
///
/// By default the controllers use the first cell next to the lower y boundary
/// (the "upstream" end of a half-domain with one target). If a measurement
/// weight is given instead, the measured value is the volume-weighted average
///
///     <f> = sum(w f dV) / sum(w dV)
///
/// over the domain interior on all processors, where w is the weight. This is
/// needed when "upstream" is not at a boundary, for example in a 1D domain with
/// a target at each end. The weight is an input expression, e.g.
///
///     exp(-(mesh:ypol / 0.02)^2)
///
/// for the cells around ypol = 0. Every processor gets the same value, so
/// controllers can update their state on all processors.
class FeedbackMeasurement {
public:
  FeedbackMeasurement() = default;

  /// Read the weight from `option`. If it is not set, usesWeight() is false
  /// and the controller keeps its original measurement location.
  FeedbackMeasurement(Options& option, const std::string& what) {
    if (!option.isSet()) {
      return;
    }
    use_weight = true;
    Mesh* mesh = bout::globals::mesh;
    const Field3D w = option
                          .doc("Weight of each cell in the measurement of " + what
                               + ": volume-weighted average where the weight is "
                                 "non-zero. If not set, the first cell next to the "
                                 "lower y boundary is used.")
                          .as<Field3D>();
    Coordinates* coord = mesh->getCoordinates();
    const auto& J = coord->J();
    const auto& dx = coord->dx();
    const auto& dy = coord->dy();
    const auto& dz = coord->dz();
    weight = 0.0;
    BoutReal local_total = 0.0;
    BOUT_FOR_SERIAL(i, w.getRegion("RGN_NOBNDRY")) {
      if (w[i] < 0.0) {
        throw BoutException("Feedback measurement weight for {:s} must be >= 0", what);
      }
      weight[i] = w[i] * J[i] * dx[i] * dy[i] * dz[i];
      local_total += weight[i];
    }
    MPI_Allreduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());
    if (!(total > 0.0)) {
      throw BoutException("Feedback measurement weight for {:s} is zero everywhere",
                          what);
    }
  }

  /// Is the measurement a weighted average (true), or the original location?
  bool usesWeight() const { return use_weight; }

  /// Volume-weighted average of f over the interior on all processors.
  /// Must be called on all processors.
  BoutReal average(const Field3D& f) const {
    BoutReal local_sum = 0.0;
    BOUT_FOR_SERIAL(i, f.getRegion("RGN_NOBNDRY")) { local_sum += weight[i] * f[i]; }
    BoutReal sum = 0.0;
    MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());
    return sum / total;
  }

private:
  bool use_weight{false};
  Field3D weight{0.0}; ///< Weight times cell volume; zero in guard cells
  BoutReal total{0.0}; ///< Sum of weight * volume over all processors
};

#endif // FEEDBACK_MEASUREMENT_H
