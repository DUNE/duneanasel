#include "duneanasel/fd/beam/Selections.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"

#include <memory>

int main(int argc, char const *argv[]) {

  TChain pch("cafmaker/cafTree");
  pch.Add(argv[1]);
  caf::StandardRecordProxy srp(&pch, "rec");

  Long64_t ents = pch.GetEntries();
  std::cout << "Input tree has " << ents << " entries." << std::endl;

  pch.GetEntry(0);

  // open a second copy of the file to copy full records out of
  TChain ch("cafmaker/cafTree");
  ch.Add(argv[1]);

  caf::StandardRecord *SR = nullptr;
  ch.SetBranchAddress("rec", &SR);

  // open the output trees
  TFile fout(argv[2], "RECREATE");
  fout.mkdir("NuMuCCLike")->cd();
  TTree *tNuMuCCLike = new TTree("cafTree", "NuMuCCLike");
  fout.mkdir("NuECCLike")->cd();
  TTree *tNuECCLike = new TTree("cafTree", "NuECCLike");
  fout.mkdir("NCLike")->cd();
  TTree *tNCLike = new TTree("cafTree", "NCLike");

  tNuMuCCLike->Branch("rec", &SR);
  tNuECCLike->Branch("rec", &SR);
  tNCLike->Branch("rec", &SR);

  for (Long64_t i = 0; i < ents; ++i) {
    pch.GetEntry(i);

    if (!srp.common.ixn.pandora.size()) {
      continue;
    }
    auto sample =
        sel::beam::fd1x8x6::numode::ApplySelection(srp.common.ixn.pandora[0]);

    switch (sample) {
    case sel::beam::kNuMuCCLike: {
      ch.GetEntry(i);
      tNuMuCCLike->Fill();
      break;
    }
    case sel::beam::kNuECCLike: {
      ch.GetEntry(i);
      tNuECCLike->Fill();
      break;
    }
    case sel::beam::kNCLike: {
      ch.GetEntry(i);
      tNCLike->Fill();
      break;
    }
    }
  }

  fout.Write();
  fout.Close();
}
