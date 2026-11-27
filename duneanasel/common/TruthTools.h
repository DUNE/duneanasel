#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/SRInteraction.h"

#include "duneanasel/common/type_traits.h"

namespace ana {

constexpr static auto kCannotFindTruth = std::numeric_limits<size_t>::max();

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline size_t GetTrueInteractionIndex(T const &nu_int) {
  if (!nu_int.truthOverlap.size()) {
    return kCannotFindTruth;
  }

  size_t idx = 0;
  double max_el = std::numeric_limits<double>::lowest();
  for (int i = 0; i < nu_int.truthOverlap.size(); ++i) {
    if (nu_int.truthOverlap[i] > max_el) {
      max_el = nu_int.truthOverlap[i];
      idx = i;
    }
  }
  if (idx >= nu_int.truth.size()) {
    return kCannotFindTruth;
  }
  return nu_int.truth[idx];
}

template <typename Tsr, typename Tixn,
          typename Csr = Proxyable_t<caf::StandardRecord, Tsr>,
          typename Cixn = Proxyable_t<caf::SRInteraction, Tixn>>
inline auto const *GetTrueInteraction(Tsr const &sr, Tixn const &ixn) {
  auto tidx = ana::GetTrueInteractionIndex(ixn);
  if (tidx == ana::kCannotFindTruth) {
    return decltype(&sr.mc.nu[tidx]){nullptr};
  }
  return &sr.mc.nu[tidx];
}

enum class TruthCat {
  kNumuCC = 0,
  kNumubarCC,
  kNueCC,
  kNuebarCC,
  kTauCC,
  kNC,
  kEntering,
  kUncategorised,
  kNoTruth
};

static std::vector<std::string> const TruthCategories = {
    "NumuCC", "NumubarCC", "NueCC",         "NuebarCC", "TauCC",
    "NC",     "Entering",  "Uncategorised", "NoTruth",
};

template <typename T, typename C = Proxyable_t<caf::SRTrueInteraction, T>>
inline TruthCat GetTrueCategory(
    T const *tixn,
    std::function<bool(caf::SRVector3D const &)> const &InDetTrueVtxVolume) {
  if (!tixn) {
    return TruthCat::kNoTruth;
  }

  if (!InDetTrueVtxVolume(tixn->vtx)) {
    return TruthCat::kEntering;
  }

  if (!tixn->iscc) {
    return TruthCat::kNC;
  }

  if (std::abs(tixn->pdg) == 12) {
    return (tixn->pdg > 0) ? TruthCat::kNueCC : TruthCat::kNuebarCC;
  }

  if (std::abs(tixn->pdg) == 14) {
    return (tixn->pdg > 0) ? TruthCat::kNumuCC : TruthCat::kNumubarCC;
  }

  if (std::abs(tixn->pdg) == 16) {
    return TruthCat::kTauCC;
  }

  return TruthCat::kUncategorised;
}

enum class CCTopo {
  kCC0Meson = 0,
  kCC1Pip,
  kCC1Pi0,
  kCC1Pim,
  kCCMultiPi,
  kCCOther,
  kNoTruth
};

static std::vector<std::string> const CCTopologies = {
    "CC0Meson",  "CC1Pip",  "CC1Pi0", "CC1Pim",
    "CCMultiPi", "CCOther", "NoTruth"};

template <typename T, typename C = Proxyable_t<caf::SRTrueInteraction, T>>
inline CCTopo GetTrueCCTopology(T const *tixn) {

  if (!tixn) {
    return CCTopo::kNoTruth;
  }

  int lepton_number = 0;
  int npip = 0, npim = 0, npi0 = 0;
  // ignore low energy gammas as they may be from nuclear de-excitation
  int ngamma = 0, ngamma_lt10mev = 0;
  int nother = 0;
  for (auto const &p : tixn-> prim) {
    switch (p.pdg) {
    case 11:
      lepton_number++;
      break;
    case -11:
      lepton_number--;
      break;
    case 13:
      lepton_number++;
      break;
    case -13:
      lepton_number--;
      break;
    case 22:
      (p.p.E < 10 ? ngamma_lt10mev : ngamma)++;
      break;
    case 211:
      npip++;
      break;
    case -211:
      npim++;
      break;
    case 111:
      npi0++;
      break;
    case 2212:
    case 2112:
      break;
    default:
      if (p.pdg < 1000000000) { // skip nuclear blobs and weird stuff
        nother++;
      }
    }
  }
  int npi = npip + npim + npi0;
  nother += ngamma;

  // weird number of leptons
  if ((tixn->pdg > 0) && (lepton_number != 1)) {
    return CCTopo::kCCOther;
  }
  if ((tixn->pdg < 0) && (lepton_number != -1)) {
    return CCTopo::kCCOther;
  }
  if (nother) {
    return CCTopo::kCCOther;
  }
  if (npi) {
    if (ngamma) {
      return CCTopo::kCCOther;
    }
    if (npi > 1) {
      return CCTopo::kCCMultiPi;
    }
    if (npip) {
      return CCTopo::kCC1Pip;
    }
    if (npim) {
      return CCTopo::kCC1Pim;
    }
    if (npi0) {
      return CCTopo::kCC1Pi0;
    }
  }

  return CCTopo::kCC0Meson;
}

} // namespace ana
