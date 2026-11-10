#include <chrono>
#include <iostream>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RNTupleDS.hxx>
#include <TFile.h>
#include <TH1D.h>

#include <Math/Vector3D.h>
#include <ROOT/RVec.hxx>

#include "duneanaobj/StandardRecord/StandardRecord.h"

using V3 = ROOT::Math::XYZVector;
template <typename T>
using RVec = ROOT::VecOps::RVec<T>;

inline RVec<double> TrackLengths(const RVec<RVec<caf::SRVector3D>> &starts,
                                const RVec<RVec<caf::SRVector3D>> &ends) {
  RVec<double> lengths;
  for (std::size_t i = 0; i < starts.size(); ++i) {
    for (std::size_t j = 0; j < starts[i].size(); ++j) {
      lengths.push_back((ends[i][j] - starts[i][j]).Mag());
    }
  }
  return lengths;
}

int main(int argc, char const *argv[]) {
  //ROOT::EnableImplicitMT();

  const std::string inputPath =
      (argc > 1) ? argv[1] : "build/example_spills.root";

  auto dataframe = ROOT::RDF::FromRNTuple("caf", inputPath);

  auto columns = dataframe.Define(
      "dlpTrackLengths", TrackLengths,
      {"rec.common.ixn.dlp.part.dlp.start", "rec.common.ixn.dlp.part.dlp.end"});
  auto dlpHistogram = columns.Histo1D<RVec<double>>(
      {"hTrackLengthDlp", "DLP track lengths;Length;Particles", 100, 0.0, 1000.0},
      "dlpTrackLengths");

  const auto loop_start = std::chrono::steady_clock::now();
  dlpHistogram.GetValue();
  auto const loop_elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - loop_start);
  std::cout << "Loop time: " << loop_elapsed.count() << " seconds"
            << std::endl;

  TFile outputFile("alltracklen_rdf.root", "RECREATE");
  dlpHistogram->Write();
  outputFile.Close();

  return 0;
}


