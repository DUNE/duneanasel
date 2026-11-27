#pragma once

#include "duneanasel/common/RecoTools.h"
#include "duneanasel/common/pdg.h"

#include "duneanasel/nd/ndlar/Geometry.h"

#include "duneanasel/nd/ndlar/numucc/Common.h"

namespace sel::beam::ndlar::numucc {

constexpr double kGeomContMargin_cm = 5.;
constexpr double kWallMargin_cm = 5.;

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool SPINEContained(T const &ixn) {
  for (auto const &part : ixn.part.dlp) {
    if (part.primary != 1) {
      continue;
    }
    if (ana::IsPrimaryWithPDG(part, {ana::pdg::kMuon, ana::pdg::kAMuon})) {
      if (part.contained != 1 &&
          !(PartcleStartContained(part, kWallMargin_cm) &&
            PartcleEscapesDownstream(part, 0))) {
        return false;
      }
    } else {
      if (part.contained != 1) {
        return false;
      }
    }
  }
  return true;
}

// selection using SPINE containment
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline Sample ApplyInclusiveSelectionV1_SPINE(T const &ixn) {
  if (!Is1MuWithVtxInFV(ixn)) {
    return kRejected;
  }
  if (!SPINEContained(ixn)) {
    return kRejected;
  }
  return kInclusive;
}

inline double DistanceToOuterWall(caf::SRVector3D const &xyz) {
  static auto const NDLArAABB = geo::ndlar::GetNDLArAABB();
  return NDLArAABB.DistanceToClosestAAWall(xyz);
}

// preliminary geomatrical containment
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool GeomContained(T const &ixn) {
  for (auto const &part : ixn.part.dlp) {
    if (part.primary != 1) {
      continue;
    }
    if (ana::IsPrimaryWithPDG(part, {ana::pdg::kMuon, ana::pdg::kAMuon})) {
      if (!(PartcleStartContained(part, kWallMargin_cm) &&
            PartcleEscapesDownstream(part, 0))) {
        return false;
      }
    } else {
      if (std::isnan(part.end.x) || std::isinf(part.end.x)) {
        return false;
      }
      if (DistanceToOuterWall(part.end) < kGeomContMargin_cm) {
        return false;
      }
    }
  }
  return true;
}

// selection using distance to wall (geometry) containment
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline Sample ApplyInclusiveSelectionV1_Geom(T const &ixn) {
  if (!Is1MuWithVtxInFV(ixn)) {
    return kRejected;
  }
  if (!GeomContained(ixn)) {
    return kRejected;
  }
  return kInclusive;
}

} // namespace sel::beam::ndlar::numucc
