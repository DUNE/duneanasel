#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleProcessor.hxx>

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include <TH1D.h>
#include "TFile.h"

#include <chrono>
#include <memory>

int main(int argc, char const *argv[]) {
  auto model = ROOT::RNTupleModel::Create();
  auto SR = model->MakeField<caf::StandardRecord>("rec");
  auto reader = ROOT::RNTupleReader::Open(std::move(model), "caf", argv[1]);

  TH1D alltracklen("h_track_lengths", "Track Lengths", 100, 0, 100);

  auto const loop_start = std::chrono::steady_clock::now();
  for (auto entryID : *reader) {
    reader->LoadEntry(entryID);

    for (auto const &nd_int : SR->common.ixn.dlp) {
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

  TFile hout("alltracklen_rntuple.root", "RECREATE");
  alltracklen.Write();
  hout.Close();

  return 0;
}