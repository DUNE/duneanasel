#include "duneanasel/nd/ndlar/Selections.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"

#include <iostream>
#include <memory>

int main(int argc, char const *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <caf file(s)>" << std::endl;
    return 1;
  }

  TChain ch("cafTree");

  for(int i = 1; i < argc; ++i){
    ch.Add(argv[i]);
  }

  caf::StandardRecord *SR = nullptr;
  if (!ch.GetBranch("rec")) {
    std::cerr << "Error: branch 'rec' not found in cafTree."
              << " Check input file content and loaded dictionaries." << std::endl;
    return 2;
  }
  ch.SetBranchAddress("rec", &SR);

  Long64_t ents = ch.GetEntries();
  if (ents <= 0) {
    std::cerr << "Error: no entries found in cafTree." << std::endl;
    return 3;
  }

  if (ch.GetEntry(0) <= 0 || SR == nullptr) {
    std::cerr << "Error: failed to read first entry from cafTree rec branch."
              << " Dictionaries may be missing at runtime." << std::endl;
    return 4;
  }

  TTree* numu_x_numu = new TTree("caf", "caf");
  TTree* numubar_x_numubar = new TTree("caf", "caf");
  TTree* nue_x_nue = new TTree("caf", "caf");
  TTree* nuebar_x_nuebar = new TTree("caf", "caf");

  double Ev;
  double Ev_reco;
  double Elep_reco;
  int mode;
  int isCC;
  double reco_q;
  int nuPDGunosc;
  int nuPDG;
  double BeRPA_A_cvwgt = 1;

  for (auto t : {numu_x_numu, numubar_x_numubar, nue_x_nue, nuebar_x_nuebar}) {
    t->Branch("Ev", &Ev, "Ev/D");
    t->Branch("Ev_reco", &Ev_reco, "Ev_reco/D");
    t->Branch("Elep_reco", &Elep_reco, "Elep_reco/D");
    t->Branch("mode", &mode, "mode/I");
    t->Branch("isCC", &isCC, "isCC/I");
    t->Branch("reco_q", &reco_q, "reco_q/D");
    t->Branch("nuPDGunosc", &nuPDGunosc, "nuPDGunosc/I");
    t->Branch("nuPDG", &nuPDG, "nuPDG/I");
    t->Branch("BeRPA_A_cvwgt", &BeRPA_A_cvwgt, "BeRPA_A_cvwgt/D");
  }

  for (Long64_t i = 0; i < ents; ++i) {
    if (ch.GetEntry(i) <= 0 || SR == nullptr) {
      continue;
    }

    for (auto const &nd_int : SR->common.ixn.dlp) {
      auto longp = ana::GetLongestParticle(nd_int);
      if (!longp || nd_int.truth.empty()) {
        continue;
      }
      
      if (sel::beam::ndlar::numode::ApplySelection(nd_int) !=
          sel::beam::ndlar::kNuMuCCLikeContained) {
        continue;
      }

      //evt.Ev = ana::GetLongestParticle(nd_int)->E;
      Ev_reco = longp->E;
      Elep_reco = longp->E;

      if (Ev_reco <= 0) {
        std::cout << "Warning: reconstructed energy for longest particle is non-positive in event " << i
                  << ". Skipping this interaction." << std::endl;
        continue;
      }

      auto truth_index = nd_int.truth[0];
      if (truth_index < 0 || static_cast<size_t>(truth_index) >= SR->mc.nu.size()) {
        continue;
      }
      auto truth_ixn = SR->mc.nu.at(truth_index);
      
      Ev = truth_ixn.E;

      mode = truth_ixn.mode;
      isCC = truth_ixn.iscc;
      nuPDGunosc = truth_ixn.pdgorig;
      nuPDG = truth_ixn.pdg;
      reco_q = truth_ixn.Q2;

      if (nuPDGunosc == 14 && nuPDG == 14) {
        numu_x_numu->Fill();
      } else if (nuPDGunosc == -14 && nuPDG == -14) {
        numubar_x_numubar->Fill();
      } else if (nuPDGunosc == 12 && nuPDG == 12) {
        nue_x_nue->Fill();
      } else if (nuPDGunosc == -12 && nuPDG == -12) {
        nuebar_x_nuebar->Fill();
      } else {
        std::cerr << "Warning: encountered unexpected neutrino PDG combination in event " << i
                  << ": nuPDGunosc = " << nuPDGunosc << ", nuPDG = " << nuPDG << std::endl;
      }
      

    }
  }

  TFile numu_x_numu_file("ND_FHC_ger_numu_x_numu_MicroProdN4p1_selected.root", "RECREATE");
  numu_x_numu->Write();
  numu_x_numu_file.Close();

  TFile numubar_x_numubar_file("ND_FHC_ger_numubar_x_numubar_MicroProdN4p1_selected.root", "RECREATE");
  numubar_x_numubar->Write();
  numubar_x_numubar_file.Close();

  TFile nue_x_nue_file("ND_FHC_ger_nue_x_nue_MicroProdN4p1_selected.root", "RECREATE");
  nue_x_nue->Write();
  nue_x_nue_file.Close();

  TFile nuebar_x_nuebar_file("ND_FHC_ger_nuebar_x_nuebar_MicroProdN4p1_selected.root", "RECREATE");
  nuebar_x_nuebar->Write();
  nuebar_x_nuebar_file.Close();

  return 0;
}
