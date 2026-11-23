#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleProcessor.hxx>

#include "duneanasel/nd/ndlar/Selections.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

int main(int argc, char const *argv[]) {
  if (argc < 3) {
    return 1;
  }

  std::string output = argv[1];

  std::vector<ROOT::Experimental::RNTupleOpenSpec> specs;
  for (const auto &f : std::vector<std::string>(argv + 2, argv + argc))
    specs.push_back({"caf", f});

  auto processor = ROOT::Experimental::RNTupleProcessor::CreateChain(specs);

  auto model = ROOT::RNTupleModel::Create();
  auto Ev_reco = model->MakeField<double>("Ev_reco");
  auto Elep_reco = model->MakeField<double>("Elep_reco");
  auto reco_q = model->MakeField<int>("reco_q");
  auto Ev = model->MakeField<double>("Ev");
  auto mode = model->MakeField<int>("mode");
  auto isCC = model->MakeField<bool>("isCC");
  auto nuPDGunosc = model->MakeField<int>("nuPDGunosc");
  auto nuPDG = model->MakeField<int>("nuPDG");
  auto BeRPA_A_cvwgt = model->MakeField<double>("BeRPA_A_cvwgt");

  auto xsec_systs = model->MakeField<std::vector<std::vector<double>>>("xsec_systs");
  auto flux_systs = model->MakeField<std::vector<std::vector<double>>>("flux_systs");

  *BeRPA_A_cvwgt = 1.0;

  auto writer = ROOT::RNTupleWriter::Recreate(std::move(model), "cafSkim", output);

  int total = 0;
  int selected = 0;
  
  for (const auto &entry: *processor) {
    const auto sr = *entry.GetPtr<caf::StandardRecord>("rec");

    for (const auto ixn: sr.common.ixn.dlp) {
      total++;
      auto classification = sel::beam::ndlar::numode::ApplySelectionV0p1(ixn);
      if (classification == sel::beam::ndlar::kRejected) {
        continue;
      }
      selected++;

      auto longp = ana::GetLongestParticle(ixn);
      if (!longp) {
        std::cerr << "No longest particle for interaction" << std::endl;
        continue;
      }

      if (ixn.truth.size() < 1) {
        std::cerr << "No truth match for interaction" << std::endl;
        continue;
      }

      const auto true_ixn = sr.mc.nu[ixn.truth.at(0)];

      // Ev_reco set to sum of all particle energies
      *Ev_reco = 0.0;
      for (auto part : ixn.part.dlp) {
        *Ev_reco += part.E;
      }

      *Elep_reco = longp->E;
      *reco_q = true_ixn.modq; // cheat
    
      *Ev = true_ixn.E;
      *mode = true_ixn.mode;
      *isCC = true_ixn.iscc;
      *nuPDGunosc = true_ixn.pdgorig;
      *nuPDG = true_ixn.pdg;

      xsec_systs->clear();
      for (const auto &syst: true_ixn.xsec_systs) {
        xsec_systs->emplace_back(syst.weights.begin(), syst.weights.end());
      }

      flux_systs->clear();
      for (const auto &syst: true_ixn.flux_systs) {
        flux_systs->emplace_back(syst.weights.begin(), syst.weights.end());
      }

      writer->Fill();
    }
  }

  std::cout << "Selected " << selected << " of " << total << " interactions" << std::endl;
}