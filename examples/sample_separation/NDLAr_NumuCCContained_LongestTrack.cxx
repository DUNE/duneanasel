#include "duneanasel/nd/ndlar/Selections.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"
#include "TH1D.h"

#include <memory>

int main(int argc, char const *argv[]) {

  TChain ch("cafTree");

  ch.Add(argv[1]);

  caf::StandardRecord *SR = nullptr;
  ch.SetBranchAddress("rec", &SR);

  Long64_t ents = ch.GetEntries();
  std::cout << "Input tree has " << ents << " entries." << std::endl;

  ch.GetEntry(0);

  TFile hout("mutracklen.root", "RECREATE");
  TH1D mutracklen("mutracklen", ";Track Length [cm]; Count", 100, 0, 1000);
  TH1D alltracklen("alltracklen", ";Track Length [cm]; Count", 100, 0, 1000);
  for (Long64_t i = 0; i < ents; ++i) {
    ch.GetEntry(i);

    std::cout << "Event: " << i << " with " << SR->common.ixn.dlp.size()
              << " dlp interactions." << std::endl;
    int j = 0;
    for (auto const &nd_int : SR->common.ixn.dlp) {
      std::cout << "  Interaction: " << j++ << " | InFV: "
                << (sel::beam::ndlar::InFV(nd_int) ? "true" : "false")
                << " | AllPrimaryParticlesContained: "
                << (sel::beam::ndlar::AllPrimaryParticlesContained(nd_int)
                        ? "true"
                        : "false");
      auto longp = ana::GetLongestParticle(nd_int);
      if (longp) {
        std::cout << " | LongestParticle: { pid: " << longp->pdg
                  << ", len: " << ana::ParticleLength(*longp) << " cm"
                  << std::endl;
      } else {
        std::cout << " | LongestParticle: null" << std::endl;
      }

      if (sel::beam::ndlar::numode::ApplySelection(nd_int) !=
          sel::beam::kNuMuCCLikeContained) {
        continue;
      } else {
        std::cout << "  => Selected!" << std::endl;
      }

      mutracklen.Fill(ana::ParticleLength(*ana::GetLongestParticle(nd_int)));

      for (auto const &p : nd_int.part.dlp) {
        alltracklen.Fill(ana::ParticleLength(p));
      }

    }
  }

  hout.Write();
}
