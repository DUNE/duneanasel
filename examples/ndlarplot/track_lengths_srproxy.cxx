#include "duneanasel/nd/ndlar/Selections.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

#include <chrono>
#include <memory>

int main(int argc, char const *argv[]) {

  // TChain ch("cafTree");
  
  // for(int i = 1; i < argc; ++i){
  //   ch.Add(argv[i]);
  // }
  
  // caf::StandardRecordProxy srp(&ch, "rec");
  
  // Long64_t ents = ch.GetEntries();
  // std::cout << "Input tree has " << ents << " entries." << std::endl;

  // ch.GetEntry(0);

  TFile fin(argv[1], "READ");
  TTree *tr = (TTree*)fin.Get("cafTree");
  caf::StandardRecordProxy srp(tr, "rec");
  Long64_t ents = tr->GetEntries();
  std::cout << "Input tree has " << ents << " entries." << std::endl;

  TH1D alltracklen("alltracklen", ";Track Length [cm]; Count", 100, 0, 1000);

  auto const loop_start = std::chrono::steady_clock::now();
  for (Long64_t i = 0; i < ents; ++i) {
    tr->LoadTree(i);
    //tr->GetEntry(i);

    int j = 0;
    for (auto const &nd_int : srp.common.ixn.dlp) {
      for (auto const &p : nd_int.part.dlp) {
        double length = std::hypot(
          p.end.x - p.start.x, 
          p.end.y - p.start.y,
          p.end.z - p.start.z
        );
        alltracklen.Fill(length);
      }
    }
  }
  auto const loop_elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - loop_start);
  std::cout << "Loop time: " << loop_elapsed.count() << " seconds"
            << std::endl;

  TFile hout("alltracklen_srproxy.root", "RECREATE");
  alltracklen.Write();
  hout.Close();
}
