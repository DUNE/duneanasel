#include "duneanasel/skimmers/MaCh3Winter26SkimModels.h"

#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "duneanasel/common/RecoTools.h"
#include "duneanasel/common/TruthTools.h"

#include "duneanasel/nd/ndlar/numucc/0pi/Selections.h"
#include "duneanasel/nd/ndlar/numucc/inc/Selections.h"

#include "duneanasel/fd/beam/Selections.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include "ROOT/RNTupleWriter.hxx"

#include <filesystem>
#include <string>
#include <vector>

struct DetSel {
  std::string name;
  ana::RecoBranch reco_stack;
  std::vector<std::string> selection_categories;
  std::function<bool(caf::SRVector3D const &)> InDetTrueVtxVolume;
  std::function<int(caf::SRInteraction const &ixn)> Selection;
};

DetSel NDLAr{"NDLAr",
             ana::RecoBranch::kSPINE,
             {"not_selected", "numucc_0pi", "numucc_inc"},
             [](caf::SRVector3D const &xyz) -> bool {
               return geo::ndlar::PointInNDLArCroppedActiveVolume(xyz);
             },
             [](caf::SRInteraction const &ixn) -> int {
               int sample_bits = 0;
               sample_bits +=
                   (sel::beam::ndlar::numucc::ApplyInclusiveSelectionV1_Geom(
                        ixn) == sel::beam::ndlar::numucc::kInclusive)
                   << 1;

               sample_bits += (sel::beam::ndlar::numucc::Apply0PiSelection(
                                   ixn) == sel::beam::ndlar::numucc::k0Pi)
                              << 2;

               return sample_bits;
             }};

DetSel FDHD{"FDHD",
            ana::RecoBranch::kPandora,
            {"not_selected", "numucc_inc", "nuecc_inc", "nc_inc"},
            [](caf::SRVector3D const &xyz) -> bool {
              return sel::beam::FDHD::InFV(xyz.x, xyz.y, xyz.z);
            },
            [](caf::SRInteraction const &ixn) -> int { return true; }};

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

  *meta.verinfo.duneanaobj_input = gl->first_version;
  *meta.verinfo.duneanaobj_process = gl->latest_version;
  *meta.verinfo.duneanasel = DUNEANASEL_VERSION;

  *meta.truth.categories = ana::TruthCategories;
  *meta.truth.cctopologies = ana::CCTopologies;
  *meta.truth.nixn_cctopo = std::vector<std::vector<int>>(
      ana::TruthCategories.size(),
      std::vector<int>(ana::CCTopologies.size(), 0));

  // currently broken in the inputs I think
  *meta.processed.swap_pot = {1, 0, 0};
  *meta.processed.nspills = 0;
  *meta.processed.nrecoixn = 0;
  *meta.processed.nrecoixn_cctopo = std::vector<std::vector<int>>(
      ana::TruthCategories.size(),
      std::vector<int>(ana::CCTopologies.size(), 0));

  for (auto const &xp : gl->wgts.xsec_params) {
    meta.syst.param_names->push_back(xp.name);
    meta.syst.param_knots->push_back(xp.vals);
    std::cout << "Expecting to find " << meta.syst.param_knots->back().size()
              << " knots for parameter: " << meta.syst.param_names->back()
              << std::endl;
  }
  for (auto const &xp : gl->wgts.flux_params) {
    meta.syst.param_names->push_back(xp.name);
    meta.syst.param_knots->push_back(xp.vals);
    std::cout << "Expecting to find " << meta.syst.param_knots->back().size()
              << " knots for parameter: " << meta.syst.param_names->back()
              << std::endl;
  }

  cafTr->GetEntry(0);
  auto &det_sel = GetDetectorSelection(*sr);
  *meta.det.name = det_sel.name;

  *meta.reco.path = ana::to_string(det_sel.reco_stack);
  *meta.reco.selection_categories = det_sel.selection_categories;

  size_t filled = 0;

  for (int i = 0; i < cafTr->GetEntries(); ++i) {
    cafTr->GetEntry(i);
    (*meta.processed.nspills)++;

    for (size_t tj = 0; tj < sr->mc.nu.size(); ++tj) {
      auto const *tixn = &sr->mc.nu.at(tj);

      auto tcat = static_cast<int>(
          ana::GetTrueCategory(tixn, det_sel.InDetTrueVtxVolume));
      auto tcctopo = static_cast<int>(ana::GetTrueCCTopology(tixn));

      meta.truth.nixn_cctopo->at(tcat)[tcctopo]++;
    }

    auto const &ixns = ana::GetInteractions(*sr, det_sel.reco_stack);

    for (size_t ij = 0; ij < ixns.size(); ++ij) {
      auto const &ixn = ixns.at(ij);
      (*meta.processed.nrecoixn)++;

      auto const *tixn = ana::GetTrueInteraction(*sr, ixn);

      *event.truth.category = static_cast<int>(
          ana::GetTrueCategory(tixn, det_sel.InDetTrueVtxVolume));
      *event.truth.cctopology = static_cast<int>(ana::GetTrueCCTopology(tixn));
      meta.processed.nrecoixn_cctopo->at(
          *event.truth.category)[*event.truth.cctopology]++;

      *event.reco.sample = det_sel.Selection(ixn);
      if (!*event.reco.sample) {
        continue;
      }

      *event.sr_entry = i;
      *event.ixn_entry = ij;
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

      *event.truth.is_cc = false;
      *event.truth.mode = 0;
      *event.truth.vtx_pos_cm = {-9999, -9999, -9999};
      *event.truth.pdgnu = 0;
      *event.truth.pdgnu_original = 0;
      *event.truth.enu_GeV = 0;
      event.truth.param_responses->clear();
      if (tixn) {
        *event.truth.is_cc = tixn->iscc;

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
                << filled << "/" << *meta.processed.nrecoixn
                << " skimmed events) spills from the input" << std::endl;
    }
  }

  meta_tuple->Fill();

  outf->Write();
}
