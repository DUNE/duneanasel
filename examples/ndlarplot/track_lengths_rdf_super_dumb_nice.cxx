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

inline RVec<double> TrackLengths(const RVec<caf::SRInteraction> &ixns) {
  RVec<double> lengths;
  for (const auto &ixn : ixns) {
    for (const auto &part : ixn.part.dlp) {
      lengths.push_back((part.end - part.start).Mag());
    }
  }
  return lengths;
}

int main(int argc, char const *argv[]) {
  const std::string inputPath =
      (argc > 1) ? argv[1] : "build/example_spills.root";

  auto dataframe = ROOT::RDF::FromRNTuple("caf", inputPath);

  auto columns = dataframe.Define("dlpTrackLengths", TrackLengths, {"rec.common.ixn.dlp"});
  auto dlpHistogram = columns.Histo1D<RVec<double>>(
      {"hTrackLengthDlp", "DLP track lengths;Length;Particles", 100, 0.0, 1000.0},
      "dlpTrackLengths");
  dlpHistogram.GetValue();

  TFile outputFile("alltracklen_rdf.root", "RECREATE");
  dlpHistogram->Write();
  outputFile.Close();

  return 0;
}


