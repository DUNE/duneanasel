#include "duneanasel/nd/ndlar/Selections.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include <ROOT/RDataFrame.hxx>

#include <ROOT/RVec.hxx>
template <typename T>
using RVec = ROOT::VecOps::RVec<T>;

inline RVec<caf::SRTrueInteraction> TruthMatchedInteraction(const RVec<caf::SRInteraction> &reco_ixns, const RVec<caf::SRTrueInteraction> &true_ixns) {
  RVec<caf::SRTrueInteraction> truth_matched_ixns;
  truth_matched_ixns.reserve(reco_ixns.size());
  
  for (const auto &reco_ixn : reco_ixns) {
    if (reco_ixn.truth.empty()) {
      truth_matched_ixns.push_back(caf::SRTrueInteraction{});
      continue;
    }
    const auto &true_ixn = true_ixns.at(reco_ixn.truth.at(0));
    truth_matched_ixns.push_back(true_ixn);
  }
  return truth_matched_ixns;
}

// Sum the energies of all particles in each interaction
inline RVec<double> RecoNeutrinoEnergy(const RVec<caf::SRInteraction> &ixns) {
  RVec<double> energies;
  energies.reserve(ixns.size());

  for (const auto &ixn : ixns) {
    double energy = 0.0;
    for (const auto &part: ixn.part.dlp)
      energy += part.E;
    energies.push_back(energy);
  }
  return energies;
}

inline RVec<double> TrueNeutrinoEnergy(const RVec<caf::SRInteraction> &reco_ixns, const RVec<caf::SRTrueInteraction> &true_ixns) {
  RVec<double> energies;
  energies.reserve(reco_ixns.size());
  
  for (const auto &reco_ixn : reco_ixns) {
    if (reco_ixn.truth.empty()) {
      energies.push_back(-1.0);
      continue;
    }
    const auto &true_ixn = true_ixns.at(reco_ixn.truth.at(0));
    energies.push_back(true_ixn.E);
  }
  return energies;
}

inline RVec<double> RecoELep(const RVec<caf::SRInteraction> &ixns) {
  RVec<double> energies;
  energies.reserve(ixns.size());

  for (const auto &ixn : ixns) {
    double energy = -1;
    for (const auto &part : ixn.part.dlp) {
      bool is_lepton = (std::abs(part.pdg) == 11 || std::abs(part.pdg) == 13);
      if (part.primary && is_lepton) {
        energy = part.E;
        break;
      }
    }
    energies.push_back(energy);
  }
  return energies;
}

inline RVec<int> Mode(const RVec<caf::SRInteraction> &reco_ixns, const RVec<caf::SRTrueInteraction> &true_ixns) {
  RVec<int> modes;
  modes.reserve(reco_ixns.size());
  
  for (const auto &reco_ixn : reco_ixns) {
    if (reco_ixn.truth.empty()) {
      modes.push_back(-1);
      continue;
    }
    const auto &true_ixn = true_ixns.at(reco_ixn.truth.at(0));
    modes.push_back(true_ixn.mode);
  }
  return modes;
}

inline RVec<int> isCC(const RVec<caf::SRInteraction> &reco_ixns, const RVec<caf::SRTrueInteraction> &true_ixns) {
  RVec<int> is_cc;
  is_cc.reserve(reco_ixns.size());
  
  for (const auto &reco_ixn : reco_ixns) {
    if (reco_ixn.truth.empty()) {
      is_cc.push_back(-1);
      continue;
    }
    const auto &true_ixn = true_ixns.at(reco_ixn.truth.at(0));
    is_cc.push_back(true_ixn.iscc);
  }
  return is_cc;
}

inline RVec<int> RecoQ(const RVec<caf::SRInteraction> &reco_ixns, const RVec<caf::SRTrueInteraction> &true_ixns) {
  RVec<int> recoq;
  recoq.reserve(reco_ixns.size());
  
  for (const auto &reco_ixn : reco_ixns) {
    if (reco_ixn.truth.empty()) {
      recoq.push_back(-1);
      continue;
    }
    const auto &true_ixn = true_ixns.at(reco_ixn.truth.at(0));
    recoq.push_back(true_ixn.modq);
  }
  return recoq;
}

inline RVec<int> nuPDGunosc(const RVec<caf::SRInteraction> &reco_ixns, const RVec<caf::SRTrueInteraction> &true_ixns) {
  RVec<int> nupdgunosc;
  nupdgunosc.reserve(reco_ixns.size());
  
  for (const auto &reco_ixn : reco_ixns) {
    if (reco_ixn.truth.empty()) {
      nupdgunosc.push_back(-1);
      continue;
    }
    const auto &true_ixn = true_ixns.at(reco_ixn.truth.at(0));
    nupdgunosc.push_back(true_ixn.pdgorig);
  }
  return nupdgunosc;
}

inline RVec<int> nuPDG(const RVec<caf::SRInteraction> &reco_ixns, const RVec<caf::SRTrueInteraction> &true_ixns) {
  RVec<int> nupdg;
  nupdg.reserve(reco_ixns.size());
  
  for (const auto &reco_ixn : reco_ixns) {
    if (reco_ixn.truth.empty()) {
      nupdg.push_back(-1);
      continue;
    }
    const auto &true_ixn = true_ixns.at(reco_ixn.truth.at(0));
    nupdg.push_back(true_ixn.pdg);
  }
  return nupdg;
}

// inline RVec<double> RecoELep(const RVec<caf::SRInteraction> &ixns) {
//   RVec<double> energies;
//   for (const auto &ixn : ixns) {
//     double energy = 0.0;
//     for (const auto &part : ixn.part.dlp) {
//       bool is_lepton = std::abs(part.pdg) == 11 || std::abs(part.pdg) == 13 || std::abs(part.pdg) == 15;
//       if (is_lepton & part.primary) {
//         energy = part.E;
//         break;
//       }
//     }
//     energies.push_back(energy);
//   }
//   return energies;
// }



int main(int argc, char const *argv[]) {

  ROOT::RDataFrame caf("caf", argv[1]);

  auto out = caf
    .Define("Ev_reco", RecoNeutrinoEnergy, {"rec.common.ixn.dlp"})
    .Define("Elep_reco", RecoELep, {"rec.common.ixn.dlp"})
    .Define("reco_q", RecoQ, {"rec.common.ixn.dlp", "rec.mc.nu"})
    .Define("Ev", TrueNeutrinoEnergy, {"rec.common.ixn.dlp", "rec.mc.nu"})
    .Define("mode", Mode, {"rec.common.ixn.dlp", "rec.mc.nu"})
    .Define("isCC", isCC, {"rec.common.ixn.dlp", "rec.mc.nu"})
    .Define("nuPDGunosc", nuPDGunosc, {"rec.common.ixn.dlp", "rec.mc.nu"})
    .Define("nuPDG", nuPDG, {"rec.common.ixn.dlp", "rec.mc.nu"})
    .Define("BeRPA_A_cvwgt", "1.0");
  
  //auto out = caf.Define("Ev", NeutrinoEnergies, {"rec.common.ixn.dlp"});
  //out = out.Define("Ev_reco", "Ev");

  ROOT::RDF::RSnapshotOptions opts;
  opts.fOutputFormat = ROOT::RDF::ESnapshotOutputFormat::kRNTuple;
  out.Snapshot("cafSkim", "tracks.root", {"Ev_reco", "Ev", "Elep_reco", "mode", "isCC", "reco_q", "nuPDGunosc", "nuPDG", "BeRPA_A_cvwgt"}, opts);
}