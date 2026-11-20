#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

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
  std::function<bool(caf::SRInteractionProxy const *)> cut;
};

struct DetSel {
  std::function<bool(double, double, double)> InDetTrueVtxVolume;
  std::vector<NamedCut> CutFlow;

  int PassesCuts(caf::SRInteractionProxy const *ixn) {
    int cut_mask = 0;
    for (size_t lvl = 0; lvl < CutFlow.size(); ++lvl) {
      cut_mask += (!CutFlow[lvl].cut(ixn) << lvl);
    }
    return cut_mask;
  }
};

static constexpr bool kSPINEReco = true;
static constexpr bool kPandoraReco = false;

template <bool isSPINE = true>
bool Has1PrimRecoMulike(caf::SRInteractionProxy const *ixn) {
  size_t nprimrecomu = 0;

  auto const *parts = &ixn->part.dlp;
  if constexpr (!isSPINE) {
    parts = &ixn->part.pandora;
  }

  for (auto const &p : *parts) {
    if ((!p.primary) || (p.pdg != 13)) {
      continue;
    }
    nprimrecomu += 1;
  }
  return nprimrecomu == 1;
}

DetSel NDLAr{[](double x, double y, double z) {
               return ana::NDLArGeo::IsInOuterFiducialVolume(x, y, z);
             },
             {
                 NamedCut{"Vertex in FV",
                          [](caf::SRInteractionProxy const *ixn) {
                            return ana::NDLArGeo::IsInCombinedNDLArFV(
                                ixn->vtx.x, ixn->vtx.y, ixn->vtx.z, 50, 50, 50,
                                50, 150, 50, 10);
                          }},
                 NamedCut{"Has 1 Primary Muon-like track",
                          Has1PrimRecoMulike<kSPINEReco>},
             }};

DetSel FDHD{
    [](double x, double y, double z) { return sel::beam::FDHD::InFV(x, y, z); },
    {
        NamedCut{"Vertex in FV",
                 [](caf::SRInteractionProxy const *ixn) {
                   return sel::beam::FDHD::InFV(*ixn);
                 }},
        NamedCut{"Has 1 Primary Muon-like track",
                 Has1PrimRecoMulike<kPandoraReco>},
    }};

caf::SRTrueInteractionProxy *
GetTrueInteraction(caf::SRInteractionProxy const *ixn) {
  auto tidx = ana::GetTrueInteractionIndex(*ixn);
  if (tidx == ana::kCannotFindTruth) {
    return nullptr;
  }
  auto srp = ixn->Ancestor<caf::StandardRecordProxy>();
  return &srp->mc.nu[tidx];
}

enum TruthCategory {
  kNumuCC = 0,
  kNC,
  kNue,
  kEntering,
  kNoTruth,
  kNCategories
};
TruthCategory GetTrueCategory(
    caf::SRTrueInteractionProxy const *tixn,
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

// these branches aren't filled
// enum Detector { nd_lar = 0, fd_hd, fd_vd };
// Detector GetDetector(caf::StandardRecordProxy &srp) {
//   if (srp.meta.nd_lar.enabled) {
//     return nd_lar;
//   }
//   if (srp.meta.fd_hd.enabled) {
//     return fd_hd;
//   }
//   if (srp.meta.fd_vd.enabled) {
//     return fd_vd;
//   }
//   throw std::runtime_error("No valid detector information found in
//   srp.meta");
// }

// DetSel &GetDetectorSelection(caf::StandardRecordProxy &srp) {
//   switch (GetDetector(srp)) {
//   case nd_lar:
//     return NDLAr;
//   case fd_hd:
//     return FDHD;
//   default: {
//     throw std::runtime_error(
//         "Can only currently process nd_lar and fd_hd files.");
//   }
//   }
// }

DetSel &GetDetectorSelection(std::string const &farg) {
  if (farg == "nd_lar") {
    return NDLAr;
  }
  if (farg == "fd_hd") {
    return FDHD;
  }

  throw std::runtime_error(
      std::string(
          "Expected the first argument to be nd_lar or fd_hd, but found ") +
      farg);
}

int main(int argc, const char *argv[]) {

  auto &det_sel = GetDetectorSelection(argv[1]);

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
  auto sr_entry = ixn_model->MakeField<uint64_t>("sr_entry");
  auto ixn_entry = ixn_model->MakeField<uint64_t>("ixn_entry");
  auto vtx_pos_cm = ixn_model->MakeField<ROOT::Math::XYZPoint>("vtx_pos_cm");
  auto enu_lep_calo_GeV = ixn_model->MakeField<double>("enu_lep_calo_GeV");
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

  std::unique_ptr<TFile> outf = std::make_unique<TFile>(argv[2], "RECREATE");

  auto meta_tuple =
      ROOT::RNTupleWriter::Append(std::move(meta_model), "meta", *outf);

  auto ixn_tuple =
      ROOT::RNTupleWriter::Append(std::move(ixn_model), "ixn", *outf);

  std::unique_ptr<TFile> fin(TFile::Open(argv[3], "READ"));
  auto tr = fin->Get<TTree>("cafTree");
  caf::StandardRecordProxy srp(tr, "rec");

  caf::SRGlobal *gl = nullptr;
  auto global = fin->Get<TTree>("globalTree");
  global->SetBranchAddress("global", &gl);
  global->GetEntry(0);
  gl->UpdateVersionInformation();

  *duneanaobj_input_version = gl->first_version;
  *duneanaobj_process_version = gl->latest_version;
  *duneanasel_version = DUNEANASEL_VERSION;

  *pot_processed = 1; // currently broken in the inputs I think
  *reco_path = "SPINE";

  for (auto const &xp : gl->wgts.xsec_params) {
    params_names->push_back(xp.name);
    params_knots->push_back(xp.vals);
  }
  for (auto const &xp : gl->wgts.flux_params) {
    params_names->push_back(xp.name);
    params_knots->push_back(xp.vals);
  }

  meta_tuple->Fill();

  size_t filled = 0;

  for (int i = 0; i < tr->GetEntries(); ++i) {
    tr->LoadTree(i);

    for (size_t j = 0; j < srp.common.ixn.dlp.size(); ++j) {
      auto const &ixn = srp.common.ixn.dlp[j];

      if (!det_sel.PassesCuts(&ixn)) {
        continue;
      }

      *sr_entry = i;
      *ixn_entry = j;
      *vtx_pos_cm = {ixn.vtx.x, ixn.vtx.y, ixn.vtx.z};
      *enu_lep_calo_GeV = ixn.Enu.lep_calo;
      *shower_e_GeV = 0;
      *sum_hadtrack_e_GeV = 0;
      for (auto const &p : ixn.part.dlp) {
        if (p.origRecoObjType == caf::kTrack) {
          if (p.primary && p.pdg) {
            *mu_e_GeV = p.E;
          } else {
            *shower_e_GeV += p.E;
          }
        } else if (p.origRecoObjType == caf::kShower) {
          *shower_e_GeV += p.E;
        }
      }

      auto const *tixn = GetTrueInteraction(&ixn);

      *true_category = GetTrueCategory(tixn, det_sel.InDetTrueVtxVolume);
      *true_mode = 0;
      *true_vtx_pos_cm = {-9999, -9999, -9999};
      *true_pdgnu = 0;
      *true_pdgnu_original = 0;
      *true_enu_GeV = 0;
      if (tixn) {
        *true_mode = tixn->mode;
        *true_vtx_pos_cm = {tixn->vtx.x, tixn->vtx.y, tixn->vtx.z};
        *true_pdgnu = tixn->pdg;
        *true_pdgnu_original = tixn->pdgorig;
        *true_enu_GeV = tixn->E;

        for (auto const &r : tixn->xsec_systs) {
          true_param_responses->emplace_back();
          for (auto const &w : r.weights) {
            true_param_responses->back().push_back(w);
          }
        }
        for (auto const &r : tixn->flux_systs) {
          true_param_responses->emplace_back();
          for (auto const &w : r.weights) {
            true_param_responses->back().push_back(w);
          }
        }
      }

      ixn_tuple->Fill();
      filled++;
    }
    if (i && !(i % 100)) {
      std::cout << "process " << i << "/" << tr->GetEntries() << " (" << filled
                << " skimmed) entries from the input" << std::endl;
    }
  }
  outf->Write();
}
