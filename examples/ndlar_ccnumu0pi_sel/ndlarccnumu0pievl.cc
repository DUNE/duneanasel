#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct NamedCut {
  std::string name;
  std::function<bool(caf::SRInteractionProxy const *)> cut;
};

struct NamedProj1D {
  std::string shortname;
  std::string xlabel;

  struct {
    int nbins;
    double low;
    double up;
  } bindef;
  std::function<double(caf::SRInteractionProxy const *)> xvar;
  int valid_from_level = 0;
};

bool kInFV(caf::SRInteractionProxy const *ixn) { return true; }

int kNShowers(caf::SRInteractionProxy const *ixn) { return 0; }

std::vector<NamedCut> CutFlow = {
    {"Vertex in FV", kInFV},
    {"NShowers == 0",
     [](caf::SRInteractionProxy const *ixn) { return kNShowers(ixn) == 0; }},
    // {"NPions == 0", kNoPions},
    // {"One Primary Muon", kOnePrimaryMuon},
    // {"Muon Track Length > 20 cm", kPrimaryMuonTrackLength},
    // {"Muon Contained in NDLAr", kPrimaryMuonContained},
    // {"Muon Leaves NDLAr at +z face", kPrimaryMuonLeavesOutBack},
    // {"All Other Tracks Contained in NDLAr", kNonMuonTracksContained},
};

std::vector<NamedProj1D> Projections1D = {
    {"vtxx",
     "Vertex X [cm]",
     {100, -100, 100}, // 10 bins from -100 to 100
     [](caf::SRInteractionProxy const *ixn) { return ixn->vtx.x; }},

    // {"mumom", "Muon Momentum [GeV/c]", {100, 0, 3}, kMuonMomentum, 3},
};

enum TruthCategory { kSignal = 0, kBkg, kNCategories };

std::string TruthCategoryName(int tc) {
  switch (TruthCategory(tc)) {
  case kSignal:
    return "Signal";
  case kBkg:
    return "Bkg";
  default:
    return "UnknownTruthCategory";
  }
}

TruthCategory kTrueCategory(caf::SRInteractionProxy const *ixn) {
  return kSignal;
}

fs::path cafdir =
    "/pnfs/dune/persistent/physicsgroups/dunendsim/abooth/nd-production/"
    "MiniProdN5/run-cafmaker/"
    "MiniProdN5p3_NDComplex_FHC.caf.full.sanddrift.spineonly/CAF/0000000";
int NCAFsToProcess = 10;

int main(int argc, const char *argv[]) {

  //[level][proj][true_category]
  std::vector<std::vector<std::vector<TH1D>>> hists1D;

  // set up all the histograms
  // make sure the outer vector is large enough
  hists1D.resize(CutFlow.size());
  for (size_t lvl = 0; lvl < CutFlow.size(); ++lvl) {
    // make sure the first inner vector is large enough
    hists1D[lvl].resize(Projections1D.size());
    for (size_t proj = 0; proj < Projections1D.size(); ++proj) {
      for (int tc = kSignal; tc < kNCategories; ++tc) {
        std::string hname = "Cutlvl" + std::to_string(lvl) + "_proj_ " +
                            Projections1D[proj].shortname + "_" +
                            TruthCategoryName(tc);
        // root histogram constructor shortcut, set the title to
        // "title;xtitle;ytitle"
        std::string htitle = "Passes level " + std::to_string(lvl) + ": " +
                             CutFlow[lvl].name + ";" +
                             Projections1D[proj].xlabel + ";count";

        // calls the TH1D constructor and puts the resulting object into the
        // hists1D vector
        hists1D[lvl][proj].emplace_back(
            hname.c_str(), htitle.c_str(), Projections1D[proj].bindef.nbins,
            Projections1D[proj].bindef.low, Projections1D[proj].bindef.up);
      }
    }
  }
  // histograms are all set up now, lets get loopy!

  // loop over files in directory
  int fi = 0;
  for (const auto &entry : fs::directory_iterator(cafdir)) {
    auto fp = entry.path();
    // if the extension isn't root, skip
    if (fp.extension().native() != ".root") {
      std::cout << "skipping path: " << fp.native() << std::endl;
      continue;
    }

    // open the file and connect the StandardRecordProxy
    std::cout << "fi: " << fi << " -> " << fp.native() << std::endl;
    std::unique_ptr<TFile> fin(TFile::Open(fp.native().c_str(), "READ"));
    auto tr = fin->Get<TTree>("cafTree");
    caf::StandardRecordProxy srp(tr, "rec");

    for (int i = 0; i < tr->GetEntries(); ++i) {
      tr->LoadTree(i);

      for (auto const &ixn : srp.common.ixn.dlp) {

        std::vector<double> projvals(Projections1D.size());
        // fill the projection cache with bad values
        std::fill_n(projvals.begin(), Projections1D.size(),
                    std::numeric_limits<double>::signaling_NaN());

        TruthCategory true_category = kTrueCategory(&ixn);

        // for a single event, loop over cuts filling relevant histograms as we
        // go
        for (size_t lvl = 0; lvl < CutFlow.size(); ++lvl) {

          // stop processing event if we don't pass a cut
          if (!CutFlow[lvl].cut(&ixn)) {
            break;
          }

          for (size_t proj = 0; proj < Projections1D.size(); ++proj) {
            // skip projections that cannot be made for this level
            if (Projections1D[proj].valid_from_level > lvl) {
              continue;
            }

            // check that we've calculated all projections for this level or
            // below
            if (Projections1D[proj].valid_from_level == lvl) {
              projvals[proj] = Projections1D[proj].xvar(&ixn);
            }

            // fill the histogram
            hists1D[lvl][proj][true_category].Fill(projvals[proj]);
          }
        }
      }
    }
    fi++;
    if (fi > NCAFsToProcess) {
      break;
    }
  }

  // write everythingout
  TFile fout("HistogramsOut.root", "RECREATE");
  for (auto &hists_cutlvls : hists1D) {
    for (auto &hists_projs : hists_cutlvls) {
      for (auto &hist : hists_projs) {
        fout.WriteTObject(&hist, hist.GetName());
      }
    }
  }
}
