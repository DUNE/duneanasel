#include "duneanasel/nd/ndlar/Selections.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

int main(int argc, char const *argv[]) {

  TChain ch("cafTree");

  for(int i = 1; i < argc; ++i){
    ch.Add(argv[i]);
  }

  caf::StandardRecord *SR = nullptr;
  ch.SetBranchAddress("rec", &SR);

  Long64_t ents = ch.GetEntries();
  std::cout << "Input tree has " << ents << " entries." << std::endl;

  ch.GetEntry(0);

  TTree trackTree("trackTree", "Track lengths");
  double trackLength = 0.0;
  trackTree.Branch("track length", &trackLength, "track length/D");

  for (Long64_t i = 0; i < ents; ++i) {
    ch.GetEntry(i);

    int j = 0;
    for (auto const &nd_int : SR->common.ixn.dlp) {
      for (auto const &p : nd_int.part.dlp) {
        double length = std::hypot(
          p.end.x - p.start.x, 
          p.end.y - p.start.y,
          p.end.z - p.start.z
        );
        trackLength = length;
        trackTree.Fill();
      }
    }
  }

  TFile fout("tracks.root", "RECREATE");
  trackTree.Write();
  fout.Close();
}