#include "duneanasel/skimmers/MaCh3Winter26SkimModels.h"

#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "duneanasel/common/RecoTools.h"
#include "duneanasel/common/TruthTools.h"
#include "duneanasel/nd/ndlar/NDLArGeometry.h"

#include "duneanasel/fd/beam/Selections.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "ROOT/RNTupleWriter.hxx"

#include <filesystem>
#include <string>
#include <vector>

struct NamedCut {
  std::string name;
  std::function<bool(caf::SRInteraction const &)> cut;
};

struct DetSel {
  ana::RecoBranch reco_stack;
  std::function<bool(double, double, double)> InDetTrueVtxVolume;
  std::vector<NamedCut> CutFlow;

  int CutBits(caf::SRInteraction const &ixn) {
    int cut_mask = 0;
    for (size_t lvl = 0; lvl < CutFlow.size(); ++lvl) {
      cut_mask += (!CutFlow[lvl].cut(ixn) << lvl);
    }
    return cut_mask;
  }
};

bool Has1PrimRecoMulike(std::vector<caf::SRRecoParticle> const &parts) {
  size_t nprimrecomu = 0;
  for (auto const &p : parts) {
    if ((!p.primary) || (p.pdg != 13)) {
      continue;
    }
    nprimrecomu += 1;
  }
  return nprimrecomu == 1;
}

constexpr static bool kSPINE = true;
constexpr static bool kPandora = true;

DetSel NDLAr{
    ana::RecoBranch::kSPINE,
    [](double x, double y, double z) {
      return ana::NDLArGeo::IsInOuterFiducialVolume(x, y, z);
    },
    {
        NamedCut{"Vertex in FV",
                 [](caf::SRInteraction const &ixn) {
                   return ana::NDLArGeo::IsInCombinedNDLArFV(
                       ixn.vtx.x, ixn.vtx.y, ixn.vtx.z, 50, 50, 50, 50, 150, 50,
                       10);
                 }},
        NamedCut{"Has 1 Primary Muon-like track",
                 [](caf::SRInteraction const &ixn) {
                   return Has1PrimRecoMulike(
                       ana::GetRecoParticles<kSPINE>(ixn));
                 }},
        NamedCut{"Is Contained",
                 [](caf::SRInteraction const &ixn) { return ixn.contained(); }},
    }};

DetSel FDHD{
    ana::RecoBranch::kPandora,
    [](double x, double y, double z) { return sel::beam::FDHD::InFV(x, y, z); },
    {NamedCut{"Vertex in FV",
              [](caf::SRInteraction const &ixn) {
                return sel::beam::FDHD::InFV(ixn.vtx.x, ixn.vtx.y, ixn.vtx.z);
              }},
     NamedCut{"Has 1 Primary Muon-like track",
              [](caf::SRInteraction const &ixn) {
                return Has1PrimRecoMulike(ana::GetRecoParticles<kPandora>(ixn));
              }},
     NamedCut{"Is Contained",
              [](caf::SRInteraction const &ixn) { return ixn.contained(); }}}};

constexpr static int kFV = 0b1;
constexpr static int kOnePrimMu = 0b10;
constexpr static int kContained = 0b100;
constexpr static int kAllCuts = 0b111;

caf::SRTrueInteraction const *
GetTrueInteraction(caf::StandardRecord const &sr,
                   caf::SRInteraction const &ixn) {
  auto tidx = ana::GetTrueInteractionIndex(ixn);
  if (tidx == ana::kCannotFindTruth) {
    return nullptr;
  }
  return &sr.mc.nu[tidx];
}

enum TruthCategory { kNumuCC = 0, kNC, kNue, kEntering, kNoTruth };

TruthCategory GetTrueCategory(
    caf::SRTrueInteraction const *tixn,
    std::function<bool(double, double, double)> const &InDetTrueVtxVolume) {
  if (!tixn) {
    return kNoTruth;
  }

  if (!InDetTrueVtxVolume(tixn->vtx.x, tixn->vtx.y, tixn->vtx.z)) {
    return kEntering;
  }

  if (!tixn->iscc) {
    return kNC;
  }

  if (!tixn->pdg) {
    return kNue;
  }

  return kNumuCC;
}

enum Detector { nd_lar = 0, fd_hd, fd_vd };
Detector GetDetector(caf::StandardRecord const &sr) {
  if (sr.meta.nd_lar.enabled || sr.meta.lar2x2.enabled) {
    return nd_lar;
  }
  if (sr.meta.fd_hd.enabled) {
    return fd_hd;
  }
  if (sr.meta.fd_vd.enabled) {
    return fd_vd;
  }
  throw std::runtime_error("No valid detector information found in sr.meta");
}

DetSel &GetDetectorSelection(caf::StandardRecord const &sr) {
  switch (GetDetector(sr)) {
  case nd_lar:
    return NDLAr;
  case fd_hd:
    return FDHD;
  default: {
    throw std::runtime_error(
        "Can only currently process nd_lar and fd_hd files.");
  }
  }
}

int main(int argc, const char *argv[]) {

  mach3winter26::Meta meta;
  mach3winter26::Event event;

  std::unique_ptr<TFile> outf = std::make_unique<TFile>(argv[1], "RECREATE");

  auto meta_tuple =
      ROOT::RNTupleWriter::Append(std::move(meta.model), "meta", *outf);
  auto ixn_tuple =
      ROOT::RNTupleWriter::Append(std::move(event.model), "event", *outf);

  std::unique_ptr<TFile> fin(TFile::Open(argv[2], "READ"));

  caf::StandardRecord *sr = nullptr;
  auto cafTr = fin->Get<TTree>("cafTree");
  cafTr->SetBranchAddress("rec", &sr);

  caf::SRGlobal *gl = nullptr;
  auto global = fin->Get<TTree>("globalTree");
  global->SetBranchAddress("global", &gl);
  global->GetEntry(0);
  gl->UpdateVersionInformation();

  *meta.duneanaobj_input_version = gl->first_version;
  *meta.duneanaobj_process_version = gl->latest_version;
  *meta.duneanasel_version = DUNEANASEL_VERSION;

  *meta.pot_processed = 1; // currently broken in the inputs I think

  for (auto const &xp : gl->wgts.xsec_params) {
    meta.param_names->push_back(xp.name);
    meta.param_knots->push_back(xp.vals);
    std::cout << "Expecting to find " << meta.param_knots->back().size()
              << " knots for parameter: " << meta.param_names->back()
              << std::endl;
  }
  for (auto const &xp : gl->wgts.flux_params) {
    meta.param_names->push_back(xp.name);
    meta.param_knots->push_back(xp.vals);
    std::cout << "Expecting to find " << meta.param_knots->back().size()
              << " knots for parameter: " << meta.param_names->back()
              << std::endl;
  }

  cafTr->GetEntry(0);
  auto &det_sel = GetDetectorSelection(*sr);
  *meta.reco_path = ana::to_string(det_sel.reco_stack);

  size_t filled = 0, allixn = 0;

  for (int i = 0; i < cafTr->GetEntries(); ++i) {
    cafTr->GetEntry(i);

    auto const &ixns = ana::GetInteractions(*sr, det_sel.reco_stack);

    for (size_t j = 0; j < ixns.size(); ++j) {
      auto const &ixn = ixns.at(j);
      allixn++;

      // bit field where failed cuts are marked as 1, so passes all == 0
      constexpr static int cut_mask = kAllCuts ^ kContained;
      if (det_sel.CutBits(ixn) & cut_mask) {
        continue;
      }

      *event.sr_entry = i;
      *event.ixn_entry = j;
      *event.reco.vtx_pos_cm = {ixn.vtx.x, ixn.vtx.y, ixn.vtx.z};
      *event.reco.shower_e_GeV = 0;
      *event.reco.sum_hadtrack_e_GeV = 0;

      for (auto const &p : ana::GetRecoParticles(ixn, det_sel.reco_stack)) {
        auto trk_shw = ana::GetTrackShower(p, det_sel.reco_stack);
        if (trk_shw == caf::kTrack) {
          if (p.primary && (p.pdg == 13)) {
            *event.reco.mu_e_GeV = p.E;
          } else {
            *event.reco.sum_hadtrack_e_GeV += p.E;
          }
        } else if (trk_shw == caf::kShower) {
          *event.reco.shower_e_GeV += p.E;
        }
      }

      *event.reco.enu_GeV =
          (*event.reco.shower_e_GeV + *event.reco.sum_hadtrack_e_GeV +
           *event.reco.mu_e_GeV);

      auto const *tixn = GetTrueInteraction(*sr, ixn);

      *event.truth.category = GetTrueCategory(tixn, det_sel.InDetTrueVtxVolume);
      *event.truth.mode = 0;
      *event.truth.vtx_pos_cm = {-9999, -9999, -9999};
      *event.truth.pdgnu = 0;
      *event.truth.pdgnu_original = 0;
      *event.truth.enu_GeV = 0;
      event.truth.param_responses->clear();
      if (tixn) {
        *event.truth.mode = tixn->mode;
        *event.truth.vtx_pos_cm = {tixn->vtx.x, tixn->vtx.y, tixn->vtx.z};
        *event.truth.pdgnu = tixn->pdg;
        *event.truth.pdgnu_original = tixn->pdgorig;
        *event.truth.enu_GeV = tixn->E;

        for (auto const &r : tixn->xsec_systs) {
          event.truth.param_responses->push_back(r.weights);
        }
        for (auto const &r : tixn->flux_systs) {
          event.truth.param_responses->push_back(r.weights);
        }
      }

      ixn_tuple->Fill();
      filled++;
    }
    if (i && !(i % 10)) {
      std::cout << "processed " << i << "/" << cafTr->GetEntries() << " ("
                << filled << "/" << allixn
                << " skimmed events) spills from the input" << std::endl;
    }
  }

  meta_tuple->Fill();

  outf->Write();
}
