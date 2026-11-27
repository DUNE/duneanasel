#pragma once

#include "duneanasel/nd/ndlar/numucc/Constants.h"

namespace sel::beam::ndlar::numucc {

enum Sample { kRejected = 0, kInclusive, k0Pi };

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool VtxInCombinedFV(T const &ixn) {
  return geo::ndlar::PointInCroppedLArVolume(
      ixn.vtx, kNDLAr_Crop_X_Hi_cm, kNDLAr_Crop_X_Lo_cm, kNDLAr_Crop_Y_Hi_cm,
      kNDLAr_Crop_Y_Lo_cm, kNDLAr_Crop_Z_Hi_cm, kNDLAr_Crop_Z_Lo_cm,
      kTPC_Crop_cm);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool Is1MuWithVtxInFV(T const &ixn) {
  if (!VtxInCombinedFV(ixn)) {
    return false;
  }
  if (ana::CountPrimariesWithPDG(ixn, {ana::pdg::kMuon, ana::pdg::kAMuon},
                                 ana::RecoBranch::kSPINE) != 1) {
    return false;
  }

  if (ana::ParticleLength_cm(*ana::GetFirstPrimaryWithPDG(
          ixn, {ana::pdg::kMuon, ana::pdg::kAMuon}, ana::RecoBranch::kSPINE)) <
      kMuon_Minimum_Track_Length_cm) {
    return false;
  }
  return true;
}

template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline bool PartcleStartContained(T const &part, double uniform_crop_cm) {

  static auto const NDLArAABB = geo::ndlar::GetNDLArAABB();
  static auto const wall_crop = geo::ndlar::MakeUniformCrop(uniform_crop_cm);

  return NDLArAABB.CroppedContains(part.start, wall_crop);
}

template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline bool PartcleEscapesDownstream(T const &p, double uniform_crop_cm) {

  static auto const NDLArAABB = geo::ndlar::GetNDLArAABB();
  // makes the end of the z box 1km downstream
  static auto const no_z_hi_crop =
      geo::ndlar::AABBCrop{uniform_crop_cm, uniform_crop_cm, uniform_crop_cm,
                           uniform_crop_cm, uniform_crop_cm, -1E5};

  return NDLArAABB.CroppedContains(p.end, no_z_hi_crop);
}

} // namespace sel::beam::ndlar::numucc
