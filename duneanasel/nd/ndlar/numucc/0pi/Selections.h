#pragma once

#include "duneanasel/common/RecoTools.h"
#include "duneanasel/common/pdg.h"

#include "duneanasel/nd/ndlar/Geometry.h"

#include "duneanasel/nd/ndlar/numucc/Common.h"

#include <limits>

namespace sel::beam::ndlar::numucc {

constexpr bool kUseSPINE = true;

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool NonPrimaryNonMuonTracksContained(T const &ixn,
                                             double Containment_Crop_cm = 5.0) {

  caf::SRVector3D min, max;
  // bounding box must include vertex
  min = ixn.vtx;
  max = ixn.vtx;

  for (const auto &part : ixn.part.dlp) {
    if ((part.primary != 1) || (std::abs(part.pdg) == 13) ||
        (ana::GetTrackShower<kUseSPINE>(part) == caf::kShower)) {
      continue;
    }

    min.x = std::min({min.x, part.start.x, part.end.x});
    max.x = std::max({max.x, part.start.x, part.end.x});

    min.y = std::min({min.y, part.start.y, part.end.y});
    max.y = std::max({max.y, part.start.y, part.end.y});

    min.z = std::min({min.z, part.start.z, part.end.z});
    max.z = std::max({max.z, part.start.z, part.end.z});
  }

  static auto const NDLArAABB = geo::ndlar::GetNDLArAABB();
  return NDLArAABB.CroppedContains(
             min, geo::ndlar::MakeUniformCrop(Containment_Crop_cm)) &&
         NDLArAABB.CroppedContains(
             max, geo::ndlar::MakeUniformCrop(Containment_Crop_cm));
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool PrimaryMuonContained(T const &ixn,
                                 double Containment_Crop_cm = 5.0) {

  auto prim_muon = ana::GetFirstPrimaryWithPDG<kUseSPINE>(
      ixn, {ana::pdg::kMuon, ana::pdg::kAMuon});

  if (!prim_muon) {
    return false;
  }

  caf::SRVector3D min, max;

  min.x = std::min(prim_muon->start.x, prim_muon->end.x);
  max.x = std::max(prim_muon->start.x, prim_muon->end.x);

  min.y = std::min(prim_muon->start.y, prim_muon->end.y);
  max.y = std::max(prim_muon->start.y, prim_muon->end.y);

  min.z = std::min(prim_muon->start.z, prim_muon->end.z);
  max.z = std::max(prim_muon->start.z, prim_muon->end.z);

  static auto const NDLArAABB = geo::ndlar::GetNDLArAABB();
  // makes the end of the z box 1km downstream
  auto const no_z_hi_crop = geo::ndlar::AABBCrop{
      Containment_Crop_cm, Containment_Crop_cm, Containment_Crop_cm,
      Containment_Crop_cm, Containment_Crop_cm, -1E5};

  return NDLArAABB.CroppedContains(
             min, geo::ndlar::MakeUniformCrop(Containment_Crop_cm)) &&
         NDLArAABB.CroppedContains(max, no_z_hi_crop);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline Sample Apply0PiSelection(T const &ixn) {
  if (!Is1MuWithVtxInFV(ixn)) {
    return kRejected;
  }
  if (ana::CountShowers<kUseSPINE>(ixn, ana::kAllParticles)) {
    return kRejected;
  }
  if (ana::CountPrimariesWithPDG<kUseSPINE>(
          ixn, {ana::pdg::kPiPlus, ana::pdg::kPiMinus})) {
    return kRejected;
  }
  if (!NonPrimaryNonMuonTracksContained(ixn)) {
    return kRejected;
  }
  if (!PrimaryMuonContained(ixn)) {
    return kRejected;
  }

  return k0Pi;
}

} // namespace sel::beam::ndlar::numucc
