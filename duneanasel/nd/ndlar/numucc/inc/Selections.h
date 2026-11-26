#pragma once

#include "duneanasel/common/RecoTools.h"
#include "duneanasel/common/type_traits.h"

#include "duneanasel/nd/ndlar/Geometry.h"

namespace sel::beam::ndlar {

enum Sample {
  kRejected = 0,
  kNuMuCCLikeContained,
  kNuECCLikeContained,
  kNuMuCCInclusive
};

// numu CC inclusive selection
namespace numuCC {

constexpr double kGeomContMargin_cm = 5.;
constexpr double kWallMargin_cm = 5.;

inline double DistToOuterWall(double x, double y, double z) {
  static auto const NDLArAABB = GetNDLArAABB();
  return NDLArAABB.DistanceToWall(x, y, z);
}

template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline bool PartcleEscapesDownstream(T const &p) {

  static auto const NDLArAABB = geo::ndlar::GetNDLArAABB();
  static auto const wall_crop = geo::ndlar::MakeUniformCrop(kWallMargin_cm);
  static auto const no_z_hi_crop = AABBCrop{0, 0, 0, 0, 0, -1E10};

  return NDLArAABB.CroppedContains(p.start, wall_crop) &&
         NDLArAABB.CroppedContains(p.end, no_z_hi_crop)
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool PassesPreSelectionV1(T const &nd_int) {
  if (!VtxInCombinedFV(nd_int)) {
    return false;
  }
  if (CountPrimariesWithPDG(nd_int, {ana::pdg::kMuon, ana::pdg::kAMuon},
                            ana::RecoBranch::SPINE) != 1) {
    return false;
  }

  auto prim_muon = ana::GetFirstPrimaryWithPDG(
      nd_int, {ana::pdg::kMuon, ana::pdg::kAMuon}, ana::RecoBranch::SPINE);

  if (ana::ParticleLength_cm(prim_muon) < kMuon_Minimum_Track_Length_cm) {
    return false;
  }
  return true;
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool SPINEContained(T const &nd_int) {
  for (auto const &p : nd_int.part.dlp) {
    if (p.primary != 1)
      continue;
    if (ana::IsPrimaryWithPDG(p, {ana::pdg::kMuon, ana::pdg::kAMuon})) {
      if (p.contained != 1 && !PartcleEscapesDownstream(p))
        return false;
    } else {
      if (p.contained != 1)
        return false;
    }
  }
  return true;
}

// selection using SPINE containment
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline sel::beam::ndlar::Sample ApplySelectionV1_SPINE(T const &nd_int) {
  if (!PassesPreSelectionV1(nd_int)) {
    return sel::beam::ndlar::kRejected;
  }
  if (!SPINEContained(nd_int)) {
    return sel::beam::ndlar::kRejected;
  }
  return sel::beam::ndlar::kNuMuCCInclusive;
}

// preliminary geomatrical containment
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool GeomContained(T const &nd_int) {
  for (auto const &p : nd_int.part.dlp) {
    if (p.primary != 1)
      continue;
    if (IsPrimaryMuon(p)) {
      if (!PartcleEscapesDownstream(p))
        return false;
    } else {
      if (std::isnan(p.end.x) || std::isinf(p.end.x))
        return false;
      if (DistToOuterWall(p.end.x, p.end.y, p.end.z) < kGeomContMargin_cm)
        return false;
    }
  }
  return true;
}

// selection using distance to wall (geometry) containment
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline sel::beam::ndlar::Sample ApplySelectionV1_Geom(T const &nd_int) {
  if (!PassesPreSelectionV1(nd_int)) {
    return sel::beam::ndlar::kRejected;
  }
  if (!GeomContained(nd_int)) {
    return sel::beam::ndlar::kRejected;
  }
  return sel::beam::ndlar::kNuMuCCInclusive;
}

} // namespace numuCC
} // namespace sel::beam::ndlar
