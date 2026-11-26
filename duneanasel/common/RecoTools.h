#pragma once

#include "duneanasel/common/type_traits.h"

#include "duneanaobj/StandardRecord/SRInteraction.h"

#include <limits>
#include <set>
#include <vector>

namespace ana {

// below is ugly, but lets us define it only once and get the enum and both
// maps https://en.wikipedia.org/wiki/X_macro
#define LIST_OF_VARIABLES                                                      \
  X(Pandora)                                                                   \
  X(SPINE)

#define X(a) k##a,
enum class RecoBranch { LIST_OF_VARIABLES };
#undef X

#define X(a)                                                                   \
  case RecoBranch::k##a:                                                       \
    return #a;
inline std::string_view const to_string(RecoBranch b) {
  switch (b) {
    LIST_OF_VARIABLES
  default: {
    throw std::runtime_error("Invalid RecoBranch specified.");
  }
  }
}
#undef X

#undef LIST_OF_VARIABLES

template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline double ParticleLength_cm(T const &p) {
  if (std::isnan(p.start.x) || std::isinf(p.start.x) || std::isnan(p.end.x) ||
      std::isinf(p.end.x)) {
    return -1.;
  }
  if (std::isnan(p.start.y) || std::isinf(p.start.y) || std::isnan(p.end.y) ||
      std::isinf(p.end.y)) {
    return -1.;
  }
  if (std::isnan(p.start.z) || std::isinf(p.start.z) || std::isnan(p.end.z) ||
      std::isinf(p.end.z)) {
    return -1.;
  }
  return std::hypot(p.end.x - p.start.x, p.end.y - p.start.y,
                    p.end.z - p.start.z);
}

template <bool isSPINE, typename T,
          typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline auto GetTrackShower(T const &p) {
  if constexpr (isSPINE) {
    switch (p.recoobj.type) {
    case caf::SRRecoBaseID::kNDLArDLPTrack:
      return caf::kTrack;
    case caf::SRRecoBaseID::kNDLArDLPShower:
      return caf::kShower;
    default:
      return caf::kUnknownRecoObj;
    }
  } else {
    switch (p.recoobj.type) {
    case caf::SRRecoBaseID::kNDLArPandoraTrack:
      return caf::kTrack;
    case caf::SRRecoBaseID::kNDLARPandoraShower:
      return caf::kShower;
    default:
      return caf::kUnknownRecoObj;
    }
  }
}

// looks silly but lets you avoid a runtime branch if you know which reco you're
// using ahead of time
template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline auto GetTrackShower(T const &p, RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE) ? GetTrackShower<true>(p)
                                      : GetTrackShower<false>(p);
}

template <bool isSPINE, typename T,
          typename C = Proxyable_t<caf::StandardRecord, T>>
inline auto const &GetInteractions(T const &sr) {
  if constexpr (isSPINE) {
    return sr.common.ixn.dlp;
  }
  return sr.common.ixn.pandora;
}

template <typename T, typename C = Proxyable_t<caf::StandardRecord, T>>
inline auto const &GetInteractions(T const &sr,
                                   RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE) ? GetInteractions<true>(sr)
                                      : GetInteractions<false>(sr);
}

template <bool isSPINE, typename T,
          typename C = Proxyable_t<caf::SRInteraction, T>>
inline auto const &GetRecoParticles(T const &ixn) {
  if constexpr (isSPINE) {
    return ixn.part.dlp;
  }
  return ixn.part.pandora;
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline auto const &GetRecoParticles(T const &ixn,
                                    RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE) ? GetRecoParticles<true>(ixn)
                                      : GetRecoParticles<false>(ixn);
}

// the auto return type here lets the compiler figure out whether to return a
//  caf::SRRecoParticle * or a caf::SRRecoParticleProxy * for a given
//  instantiation
template <bool isSPINE, typename T,
          typename C = Proxyable_t<caf::SRInteraction, T>>
inline auto const *GetLongestParticle(T const &ixn) {

  auto const &parts = GetRecoParticles<isSPINE>(ixn);

  double longest = std::numeric_limits<double>::lowest();
  size_t longest_idx = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < parts.size(); ++i) {
    const auto &p = parts.at(i);
    if (p.primary != 1) {
      continue;
    }
    double const &plen = ParticleLength(p);
    if (plen > longest) {
      longest = plen;
      longest_idx = i;
    }
  }
  return (longest_idx == std::numeric_limits<size_t>::max())
             ? nullptr
             : &parts.at(longest_idx);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline auto const *GetLongestParticle(T const &ixn,
                                      RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE) ? GetLongestParticle<true>(ixn)
                                      : GetLongestParticle<false>(ixn);
}

template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
inline bool IsPrimaryWithPDG(T const &p, std::set<int> const &pdgs) {
  return (p.primary == 1) && pdgs.count(p.pdg);
}

template <bool isSPINE, typename T,
          typename C = Proxyable_t<caf::SRInteraction, T>>
inline int CountPrimariesWithPDG(T const &nd_int, std::set<int> const &pdgs) {

  auto const &parts = GetRecoParticles<isSPINE>(ixn);

  int n = 0;
  for (auto const &p : parts) {
    n += IsPrimaryWithPDG(p, pdgs);
  }
  return n;
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline int CountPrimariesWithPDG(T const &nd_int, std::set<int> const &pdgs,
                                 RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE) ? CountPrimariesWithPDG<true>(ixn, pdgs)
                                      : CountPrimariesWithPDG<false>(ixn, pdgs);
}

/// single pdg override
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline int CountPrimariesWithPDG(T const &nd_int, int pdg,
                                 RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE)
             ? CountPrimariesWithPDG<true>(ixn, {pdg})
             : CountPrimariesWithPDG<false>(ixn, {pdg});
}

template <bool isSPINE, typename T,
          typename C = Proxyable_t<caf::SRInteraction, T>>
inline auto GetFirstPrimaryWithPDG(T const &nd_int, std::set<int> const &pdgs) {
  auto const &parts = GetRecoParticles<isSPINE>(ixn);

  std::add_pointer_t<parts::value_type> p = nullptr;

  for (auto const &p : parts) {
    if (IsPrimaryWithPDG(p, pdgs)) {
      return &p;
    }
  }
  return p;
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline auto GetFirstPrimaryWithPDG(T const &nd_int, std::set<int> const &pdgs,
                                   RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE)
             ? GetFirstPrimaryWithPDG<true>(ixn, pdgs)
             : GetFirstPrimaryWithPDG<false>(ixn, pdgs);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline auto GetFirstPrimaryWithPDG(T const &nd_int, std::set<int> pdg,
                                   RecoBranch reco = RecoBranch::kSPINE) {
  return (reco == RecoBranch::kSPINE)
             ? GetFirstPrimaryWithPDG<true>(ixn, {pdg})
             : GetFirstPrimaryWithPDG<false>(ixn, {pdg});
}

} // namespace ana
