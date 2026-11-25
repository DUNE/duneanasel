#include "duneanasel/nd/ndlar/numucc0pi/Cuts.h"
#include "duneanasel/common/TruthTools.h"
#include "duneanasel/nd/ndlar/numucc0pi/TruthTools.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "ROOT/RNTuple.hxx"
#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "Math/Point3D.h"

#include <iostream>
#include <memory>
#include <sstream>
#include <filesystem>
#include <limits>
#include <algorithm>
#include <string>
#include <cmath>

namespace fs = std::filesystem;

using namespace ana::numucc0pi;

std::string PDGToString(int pdg){
  if (pdg == 13){return "Muon";}
  else if (pdg == 11) {return "Electron";}
  else if (pdg == 211) {return "Pion";}
  else if (pdg == 2212) {return "Proton";}
  else if (pdg == 111) {return "PionZero";}
  else if (pdg == -211) {return "PionMinus";}
  else if (pdg == 22) {return "Photon";}
  else if (pdg == 310) {return "KaonZeroShort";}
  else if (pdg == 311) {return "KaonZeroLong";}
  else if (pdg == 321) {return "Kaon";}
  else if (pdg == 130) {return "KaonZero";}
  return "Particle";
}


  
std::string PrintParticle(caf::SRRecoParticleProxy const *part){
  std::stringstream ss;
  ss << "{ Reco PDG:" << PDGToString(part->pdg)
      << " Start point:  [" << part->start.x << ", " << part->start.y << ", " << part->start.z
      << "], End point:  [" << part->end.x << ", " << part->end.y << ", " << part->end.z << "] }";
  return ss.str();
}


int main(int argc, char* argv[]){

  int NMaxFilesToRead = 10;

  // --- RNTUPLE INITIALIZATION ---
  auto model = ROOT::RNTupleModel::Create();

  // Interaction-level scalar fields
  auto fEventId          = model->MakeField<int>("event_id");
  auto fInteractionId    = model->MakeField<int>("interaction_id");
  auto fVtxX             = model->MakeField<float>("vtx_x");
  auto fVtxY             = model->MakeField<float>("vtx_y");
  auto fVtxZ             = model->MakeField<float>("vtx_z");

  // Non-muon track bounding box fields
  auto fBoxMin           = model->MakeField<ROOT::Math::XYZPoint>("box_min");
  auto fBoxMax           = model->MakeField<ROOT::Math::XYZPoint>("box_max");
  auto fMuonStart        = model->MakeField<ROOT::Math::XYZPoint>("muon_start");
  auto fMuonEnd          = model->MakeField<ROOT::Math::XYZPoint>("muon_end");
  auto fMuonMom          = model->MakeField<float>("muon_mom");
  auto fPassesFixed      = model->MakeField<bool>("passes_fixed");

  // Truth category classification field (0 to 9)
  auto fTruthCategory    = model->MakeField<int>("truth_category");
  // total true signal counted from true events
  

  //true vertex information
  auto fTrueVtxX = model->MakeField<float>("true_vtx_x");
  auto fTrueVtxY = model->MakeField<float>("true_vtx_y");
  auto fTrueVtxZ = model->MakeField<float>("true_vtx_z");

  

  // cut flags
  auto fPassesNoShowers = model->MakeField<bool>("passes_no_showers");
  auto fPassesOneMuon   = model->MakeField<bool>("passes_one_muon");
  auto fPassesZeroPions = model->MakeField<bool>("passes_zero_pions");

  auto fPassesContainment = model->MakeField<bool>("passes_containment");

//variable parameter flags (not used right now)
  auto fVaryFiducial     = model->MakeField<bool>("vary_fiducial");
  auto fVaryMuLength     = model->MakeField<bool>("vary_mu_length");
  auto fVaryMuContain    = model->MakeField<bool>("vary_mu_contain");
  auto fVaryNonMuTrack   = model->MakeField<bool>("vary_nonmu_track");

  auto writer = ROOT::RNTupleWriter::Recreate(std::move(model), "myNTuple", "output_optimization_ntuple_3.root");

  //counters for particles and unknown objects
  long total_particles = 0;
  long unkown_particles = 0;
  long secondary_particles = 0;
  long long total_true_signal_in_fv = 0; //counter for efficiency denominator (true signal total)

  //fs::path cafdir = "/pnfs/dune/persistent/physicsgroups/dunendsim/abooth/nd-production/MiniProdN5/run-cafmaker/MiniProdN5p3_NDComplex_FHC.caf.full.sanddrift.spineonly/CAF/0000000";
    
  fs::path cafdir ="/exp/dune/data/users/picker24/NDMiniprod5N3/flat";

  int fi = 0;
  for (const auto& entry : fs::directory_iterator(cafdir)) {
    auto fp = entry.path();
    if(fp.extension().native() != ".root"){
      std::cout << "skipping path: " << fp.native() << std::endl;
      continue;
    }

    std::cout << "fi: " << fi << " -> " << fp.native() << std::endl;

    std::unique_ptr<TFile> fin(TFile::Open(fp.native().c_str(),"READ"));
    auto tr = fin->Get<TTree>("cafTree");
    caf::StandardRecordProxy srp(tr, "rec");

   for(int i = 0; i < tr->GetEntries(); ++i){
     tr->LoadTree(i);


     // true interaction loop for counting true signal total
    for (const caf::SRTrueInteractionProxy& mc_nu : srp.mc.nu) {
        const int true_cat = static_cast<int>(ana::numucc0pi::TruthInteractionChecker(&mc_nu));
        if (true_cat == static_cast<int>(ana::numucc0pi::kNuMuCC0PiFV_ArgTarg)) {
          total_true_signal_in_fv++;
        }
      }

      std::cout << "Event: " << i << " with " << srp.common.ixn.dlp.size()
                << " dlp interactions." << std::endl;

      //skip interactions without any particles
      if (!srp.common.ixn.dlp.size()) {
        continue;
      }


      int j = 0;
      for(const caf::SRInteractionProxy &ixn : srp.common.ixn.dlp){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 1. FIXED PRE-SELECTION CUTS (Hard Filters)
        // cuts that do NOT have variable parameters
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


       // bool passes_fixed = NoShowersInteractions(&ixn);
        //bool passes_fixed =numucc0pi::ExactlyOnePrimaryMuon(&ixn);
       // passes_fixed = passes_fixed && ana::numucc0pi::ExactlyOnePrimaryMuon(&ixn);
        //passes_fixed = passes_fixed && ana::numucc0pi::ProducesZeroPions(&ixn);

        //*fPassesFixed = passes_fixed;

        // 2. DATA struct and vertex extraction
        *fEventId          = i;
        *fInteractionId    = j;
        *fVtxX             = ixn.vtx.x;
        *fVtxY             = ixn.vtx.y;
        *fVtxZ             = ixn.vtx.z;

        *fTruthCategory = static_cast<int>(ana::numucc0pi::GetTruthInteractionType(&ixn));
        // If it fails the fixed topology, skip this interaction entirely

        //define truth interaction that has been matched to the reco intercation
        const caf::SRTrueInteractionProxy* truthInt = ana::numucc0pi::GetTruthMatch(&ixn);
        if (truthInt) {
          //extract true interaction vtx coords 
            *fTrueVtxX = truthInt->vtx.x;
            *fTrueVtxY = truthInt->vtx.y;
            *fTrueVtxZ = truthInt->vtx.z;
        } else {
            // Sentinel value for kNoMatch interactions so they don't pile up in the underflow edge bin
            *fTrueVtxX = -9999.0f;
            *fTrueVtxY = -9999.0f;
            *fTrueVtxZ = -9999.0f;
        }

        //define boolean for whether event passes any of the fixed cuts 
        bool passes_no_showers = NoShowersInteractions(&ixn);
        bool passes_one_muon   = ana::numucc0pi::ExactlyOnePrimaryMuon(&ixn);
        bool passes_zero_pions = ana::numucc0pi::ProducesZeroPions(&ixn);

        bool passes_containment = ana::numucc0pi::Nu_muCCContained(&ixn);

        //assign the boolean returns to the fields in the RNTuple
        *fPassesNoShowers = passes_no_showers;
        *fPassesOneMuon   = passes_one_muon;
        *fPassesZeroPions = passes_zero_pions;

        *fPassesContainment = passes_containment;

        bool passes_fixed = passes_no_showers && passes_one_muon && passes_zero_pions;
        *fPassesFixed = passes_fixed;
        /*
        if (!passes_fixed) {
            j++;
            writer->Fill();
            continue;
        } 
        */

        //*fVaryFiducial     = kNDLArFiducialVolume(&ixn);
        //*fVaryMuLength     = kNDLArMuTrackLength(&ixn);
        //*fVaryMuContain    = kMuonsContained(&ixn);
        //*fVaryNonMuTrack   = kNonMuonsTRACKContained(&ixn);


        // Initialize bounding box trackers for this interaction
        double xmin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::lowest();
        double ymin = std::numeric_limits<double>::max();
        double ymax = std::numeric_limits<double>::lowest();
        double zmin = std::numeric_limits<double>::max();
        double zmax = std::numeric_limits<double>::lowest();
        int non_muon_track_count = 0;

        ROOT::Math::XYZPoint muon_start_pt(0, 0, 0);
        ROOT::Math::XYZPoint muon_end_pt(0, 0, 0);
        float muon_mom_val = -1.0f;

        //3. EVENT / PARTICLE EXTRACTION LOOP
        for (const auto& part : ixn.part.dlp){

            total_particles++;

          //use new shower definiton 
            bool is_shower = (part.recoobj.type == caf::SRRecoBaseID::SRRecoBaseCollectionType:: kNDLArDLPShower);
            



            
            bool is_muon = (std::abs(part.pdg) == 13);
            //bool is_shower = (std::abs(part.pdg) == 11 || std::abs(part.pdg) == 22);
            bool is_track = !is_shower; // If it's not a shower, treat as a track
            bool is_not_muon = !is_muon;

            // Capture muon track start and end points (typically the primary muon)
            if (is_track && is_muon) {
                muon_start_pt = ROOT::Math::XYZPoint(part.start.x, part.start.y, part.start.z);
                muon_end_pt   = ROOT::Math::XYZPoint(part.end.x, part.end.y, part.end.z);
                muon_mom_val  = std::sqrt(part.p.x * part.p.x +   // <-- 4. COMPUTE MOMENTUM
                                          part.p.y * part.p.y +
                                          part.p.z * part.p.z);
            }


            if (is_track && is_not_muon) {
              non_muon_track_count++;

              // Update bounding box extremes with both start and end points
              xmin = std::min({xmin, static_cast<double>(part.start.x), static_cast<double>(part.end.x)});
              xmax = std::max({xmax, static_cast<double>(part.start.x), static_cast<double>(part.end.x)});

              ymin = std::min({ymin, static_cast<double>(part.start.y), static_cast<double>(part.end.y)});
              ymax = std::max({ymax, static_cast<double>(part.start.y), static_cast<double>(part.end.y)});

              zmin = std::min({zmin, static_cast<double>(part.start.z), static_cast<double>(part.end.z)});
              zmax = std::max({zmax, static_cast<double>(part.start.z), static_cast<double>(part.end.z)});




             // Push the populated particle into the interaction's vector

              /*
              std::cout << "Event: " << i << " ,Interaction:" << j
                        << " particle PDG: "
                        << " Contained: " << part.contained
                        << " , Primary? (1=true): "  << std::endl;
                        */
            }
        }


        // Fallback to vertex position if no non-muon tracks exist in the interaction
        if (non_muon_track_count == 0) {
            xmin = xmax = ixn.vtx.x;
            ymin = ymax = ixn.vtx.y;
            zmin = zmax = ixn.vtx.z;
        }

        // Assign bounding box corners using XYZPoint
        *fBoxMin = ROOT::Math::XYZPoint(xmin, ymin, zmin);
        *fBoxMax = ROOT::Math::XYZPoint(xmax, ymax, zmax);
        *fMuonStart  = muon_start_pt;
        *fMuonEnd    = muon_end_pt;
        *fMuonMom   = muon_mom_val;
            // --- Your coordinate boundary checks & RNTuple filling code here ---
            // Fill your RNTuple branch variables here, including:
            // - Interaction-level flags: var_fiducial, var_mu_length, etc.
            // - Particle-level kinematics: part.end.x, part.pdg, etc.

          // 4. WRITE TO RNTUPLE (Writes one entry per interaction containing all its particles)

        writer->Fill();
        j++;
      }
    }
    // break after reading NMaxFilesToRead
    //if(++fi >= NMaxFilesToRead){
    ///  break;
    
    fi++;
  }

  // Write summary count after all files and events are processed
  auto meta_model = ROOT::RNTupleModel::Create();
  auto fTotalTrueSignalFV = meta_model->MakeField<long long>("total_true_signal_in_fv");
  auto meta_writer = ROOT::RNTupleWriter::Recreate(std::move(meta_model), "meta", "output_optimization_summary_3.root");
  *fTotalTrueSignalFV = total_true_signal_in_fv;
  meta_writer->Fill();


  std::cout << "Total True signal counted: " << total_true_signal_in_fv << std::endl;
  std::cout << "SUCCESS: All data successfully written to:" << std::endl;
  std::cout << "output_optimization_ntuple_3.root" << std::endl;
}


