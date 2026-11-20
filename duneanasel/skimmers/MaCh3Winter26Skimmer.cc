#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "duneanasel/common/RecoTools.h"
#include "duneanasel/common/TruthTools.h"
#include "duneanasel/nd/ndlar/NDLArGeometry.h"

#include "duneanasel/fd/beam/Selections.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#include "ROOT/RNTupleModel.hxx"
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

  auto meta_model = ROOT::RNTupleModel::Create();
  auto duneanaobj_input_version =
      meta_model->MakeField<std::string>("duneanaobj_input_version");
  auto duneanaobj_process_version =
      meta_model->MakeField<std::string>("duneanaobj_process_version");
  auto duneanasel_version =
      meta_model->MakeField<std::string>("duneanasel_version");
  auto pot_processed = meta_model->MakeField<double>("pot_processed");
  auto reco_path = meta_model->MakeField<std::string>("reco_path");
  auto params_names =
      meta_model->MakeField<std::vector<std::string>>("param_names");
  auto params_knots =
      meta_model->MakeField<std::vector<std::vector<float>>>("param_knots");

  auto ixn_model = ROOT::RNTupleModel::Create();
  auto sr_entry = ixn_model->MakeField<int>("sr_entry");
  auto ixn_entry = ixn_model->MakeField<int>("ixn_entry");
  auto vtx_pos_cm = ixn_model->MakeField<ROOT::Math::XYZPoint>("vtx_pos_cm");
  auto enu_rec_GeV = ixn_model->MakeField<double>("enu_rec_GeV");
  auto mu_e_GeV = ixn_model->MakeField<double>("mu_e_GeV");
  auto shower_e_GeV = ixn_model->MakeField<double>("shower_e_GeV");
  auto sum_hadtrack_e_GeV = ixn_model->MakeField<double>("sum_hadtrack_e_GeV");

  auto true_category = ixn_model->MakeField<int>("true_category");
  auto true_mode = ixn_model->MakeField<int>("true_mode");
  auto true_vtx_pos_cm =
      ixn_model->MakeField<ROOT::Math::XYZPoint>("true_vtx_pos_cm");
  auto true_pdgnu = ixn_model->MakeField<double>("true_pdgnu");
  auto true_pdgnu_original =
      ixn_model->MakeField<double>("true_pdgnu_original");
  auto true_enu_GeV = ixn_model->MakeField<double>("true_enu_GeV");
  auto true_param_responses =
      ixn_model->MakeField<std::vector<std::vector<float>>>(
          "true_param_responses");

  std::unique_ptr<TFile> outf = std::make_unique<TFile>(argv[1], "RECREATE");

  auto meta_tuple =
      ROOT::RNTupleWriter::Append(std::move(meta_model), "meta", *outf);

  auto ixn_tuple =
      ROOT::RNTupleWriter::Append(std::move(ixn_model), "ixn", *outf);

  std::unique_ptr<TFile> fin(TFile::Open(argv[2], "READ"));

  caf::StandardRecord *sr = nullptr;
  auto cafTr = fin->Get<TTree>("cafTree");
  cafTr->SetBranchAddress("rec", &sr);

  caf::SRGlobal *gl = nullptr;
  auto global = fin->Get<TTree>("globalTree");
  global->SetBranchAddress("global", &gl);
  global->GetEntry(0);
  gl->UpdateVersionInformation();

  *duneanaobj_input_version = gl->first_version;
  *duneanaobj_process_version = gl->latest_version;
  *duneanasel_version = DUNEANASEL_VERSION;

  *pot_processed = 1; // currently broken in the inputs I think

  for (auto const &xp : gl->wgts.xsec_params) {
    params_names->push_back(xp.name);
    params_knots->push_back(xp.vals);
    std::cout << "Expecting to find " << params_knots->back().size()
              << " knots for parameter: " << params_names->back() << std::endl;
  }
  for (auto const &xp : gl->wgts.flux_params) {
    params_names->push_back(xp.name);
    params_knots->push_back(xp.vals);
    std::cout << "Expecting to find " << params_knots->back().size()
              << " knots for parameter: " << params_names->back() << std::endl;
  }

  cafTr->GetEntry(0);
  auto &det_sel = GetDetectorSelection(*sr);
  *reco_path = ana::to_string(det_sel.reco_stack);

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

      *sr_entry = i;
      *ixn_entry = j;
      *vtx_pos_cm = {ixn.vtx.x, ixn.vtx.y, ixn.vtx.z};
      *shower_e_GeV = 0;
      *sum_hadtrack_e_GeV = 0;

      for (auto const &p : ana::GetRecoParticles(ixn, det_sel.reco_stack)) {
        auto trk_shw = ana::GetTrackShower(p, det_sel.reco_stack);
        if (trk_shw == caf::kTrack) {
          if (p.primary && (p.pdg == 13)) {
            *mu_e_GeV = p.E;
          } else {
            *sum_hadtrack_e_GeV += p.E;
          }
        } else if (trk_shw == caf::kShower) {
          *shower_e_GeV += p.E;
        }
      }

      *enu_rec_GeV = (*shower_e_GeV + *sum_hadtrack_e_GeV + *mu_e_GeV);

      auto const *tixn = GetTrueInteraction(*sr, ixn);

      *true_category = GetTrueCategory(tixn, det_sel.InDetTrueVtxVolume);
      *true_mode = 0;
      *true_vtx_pos_cm = {-9999, -9999, -9999};
      *true_pdgnu = 0;
      *true_pdgnu_original = 0;
      *true_enu_GeV = 0;
      true_param_responses->clear();
      if (tixn) {
        *true_mode = tixn->mode;
        *true_vtx_pos_cm = {tixn->vtx.x, tixn->vtx.y, tixn->vtx.z};
        *true_pdgnu = tixn->pdg;
        *true_pdgnu_original = tixn->pdgorig;
        *true_enu_GeV = tixn->E;

        for (auto const &r : tixn->xsec_systs) {
          true_param_responses->push_back(r.weights);
        }
        for (auto const &r : tixn->flux_systs) {
          true_param_responses->push_back(r.weights);
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
