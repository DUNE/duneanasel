//R__ADD_INCLUDE_PATH(/exp/dune/app/users/fvosper/cafana-tutorialSep24/)
//R__ADD_INCLUDE_PATH(/exp/dune/app/users/fvosper/cafana-tutorialSep24/duneanasel)



//#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
//

//#include "duneanaobj/StandardRecord/StandardRecord.h"
/*
#include <iostream>
#include <memory>
#include <sstream>
#include <filesystem>
#include <limits>
#include <algorithm>



#include "duneanasel/nd/ndlar/Nu_muCC0Pi.h"
#include "duneanasel/common/TruthTools.h"
//include "duneanasel/nd/ndlar/Nu_muCC0Pi_Cuts.h"
//#include "duneanasel/lblpwgtools/CAFAna/Nu_muCC0Pi_Selection/Nu_muCC0Pi_Cuts.h"

#include "TFile.h"
#include "TTree.h"



#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include "Math/Point3D.h"

#include <memory>
#include <sstream>
#include <filesystem>

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
*/

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>


#include <iostream>
#include <memory>
#include <sstream>
#include <filesystem>
#include <limits>
#include <algorithm>
#include <string>

// 1. ROOT & RNTuple Headers
#include "TFile.h"
#include "TTree.h"
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include "Math/Point3D.h"


#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

// 2. StandardRecord Proxy Headers (MUST COME BEFORE duneanasel headers)
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

// 3. duneanasel Headers
#include "duneanasel/nd/ndlar/Nu_muCC0Pi.h"
#include "duneanasel/common/TruthTools.h"


namespace fs = std::filesystem;

using namespace ana;

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

// std::string PrintRecoObjType(int i){
//   switch(i){
//     case -1:
//       return "kUnknownRecoObj";
//     case 1:
//       return "KTrack";
//     case 2:
//       return "kShower";
//     case 3:
//       return "kHits";
//     default:
//       return "kUnknownRecoObj";
//   }
// }


//want function that takes SRParticleProxy * and returns true if track and false otherwise
//want cut that takes interacxtion adn loops ovcer particles 
//call 2nd function for cut 
/*
std::string PrintRecoObjType(caf::SRRecoBaseID::SRRecoBaseCollectionType t){
  switch(t){
    case caf::SRRecoBaseID::kNDLArDLPShower:
      return "kNDLArDLPShower";
    case caf::SRRecoBaseID::kNDLArDLPTrack:
      return "kNDLArDLPTrack";
    default:
      return std::string("Unknown:") + std::to_string(int(t));
  }
}
  */
/**
std::string PrintRecoObjType(int type_val) {
  if (type_val == caf::SRRecoBaseID::kNDLArDLPShower) return "kNDLArDLPShower";
  if (type_val == caf::SRRecoBaseID::kNDLArDLPTrack)  return "kNDLArDLPTrack";
  return std::string("Unknown:") + std::to_string(type_val);
}
*/

std::string PrintRecoObjType(int type_val) {
  if (type_val == 1) return "kNDLArDLPShower";
  if (type_val == 2) return "kNDLArDLPTrack";
  return "Unknown:" + std::to_string(type_val);
}
std::string PrintParticle(caf::SRRecoParticleProxy const *part){
  std::stringstream ss;
  
 // ss << "{ Reco PDG:" << PDGToString(part->pdg) << ", Reco Type: " << PrintRecoObjType(part->recoobj.type.GetValue())
  ss << "{ Reco PDG:" << PDGToString(part->pdg) //<< ", Reco Type: " << PrintRecoObjType(part->type)
      << " Start point:  [" << part->start.x << ", " << part->start.y << ", " << part->start.z 
      << "], End point:  [" << part->end.x << ", " << part->end.y << ", " << part->end.z << "] }";
  return ss.str();
}




//void Event_Loop(){
int main(int argc, char* argv[]){
 
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
  auto fPassesFixed      = model->MakeField<bool>("passes_fixed");

// Truth category classification field (0 to 9)
  auto fTruthCategory    = model->MakeField<int>("truth_category");

  // Variable cut optimization flags
  auto fVaryFiducial     = model->MakeField<bool>("vary_fiducial");
  auto fVaryMuLength     = model->MakeField<bool>("vary_mu_length");
  auto fVaryMuContain    = model->MakeField<bool>("vary_mu_contain");
  auto fVaryNonMuTrack   = model->MakeField<bool>("vary_nonmu_track");
 
  auto writer = ROOT::RNTupleWriter::Recreate(std::move(model), "myNTuple", "output_optimization_ntuple_2.root");


  //counters for particles and unknown objects 
  long long total_particles = 0;
  long long unkown_particles = 0;
  long long secondary_particles=0;

  fs::path cafdir = "/pnfs/dune/persistent/physicsgroups/dunendsim/abooth/nd-production/MiniProdN5/run-cafmaker/MiniProdN5p3_NDComplex_FHC.caf.full.sanddrift.spineonly/CAF/0000000";

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

      std::cout << "Event: " << i << " with " << srp.common.ixn.dlp.size()
                << " dlp interactions." << std::endl;

      //skip interactions without any particles
      if (!srp.common.ixn.dlp.size()) {
        continue;
      }

      
      int j = 0;
      for(const caf::SRInteractionProxy &ixn : srp.common.ixn.dlp){
        
        // 1. FIXED PRE-SELECTION CUTS (Hard Filters)
        // cuts that do NOT have variable parameters

        bool passes_fixed = NDLArFunctions::NoShowersInteractions(&ixn);
        //bool passes_fixed =NDLArFunctions::ExactlyOnePrimaryMuon(&ixn);
        passes_fixed = passes_fixed && NDLArFunctions::ExactlyOnePrimaryMuon(&ixn);
        passes_fixed = passes_fixed && NDLArFunctions::ProducesZeroPions(&ixn);

        *fPassesFixed = passes_fixed;

        // 2. DATA struct and vertex extraction
        *fEventId          = i;
        *fInteractionId    = j;
        *fVtxX             = ixn.vtx.x;
        *fVtxY             = ixn.vtx.y;
        *fVtxZ             = ixn.vtx.z;
        
        //*fTruthCategory    = static_cast<int>(ana::GetTruthInteractionType(&ixn));
        *fTruthCategory = static_cast<int>(ana::GetTruthInteractionType(&ixn));
        // If it fails the fixed topology, skip this interaction entirely
        if (!passes_fixed) {
          j++;
          writer->Fill ();         
          continue;
        }
      
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

        //3. EVENT / PARTICLE EXTRACTION LOOP
        for (const auto& part : ixn.part.dlp){
            
            total_particles++;

          /*
            std::string reco_type_str = PrintRecoObjType(part.recoobj.type.GetValue());
            //std::string reco_type_str = PrintRecoObjType(part.type.GetValue());
            if (reco_type_str.find("Unknown") != std::string::npos) {
              unkown_particles++;
              if (part.primary == 0){
                secondary_particles++;
              }
            }

            // Check if particle is a track
            bool is_track = (reco_type_str.find("kNDLArDLPTrack") != std::string::npos);
            bool is_muon = (std::abs(part.pdg) == 13);
            bool is_not_muon = !is_muon;
            */
            // Use PDG code to classify the particle (11=Electron, 22=Photon are showers)
            bool is_muon = (std::abs(part.pdg) == 13);
            bool is_shower = (std::abs(part.pdg) == 11 || std::abs(part.pdg) == 22);
            bool is_track = !is_shower; // If it's not a shower, treat as a track
            bool is_not_muon = !is_muon;

            // Capture muon track start and end points (typically the primary muon)
            if (is_track && is_muon) {
                muon_start_pt = ROOT::Math::XYZPoint(part.start.x, part.start.y, part.start.z);
                muon_end_pt   = ROOT::Math::XYZPoint(part.end.x, part.end.y, part.end.z);
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
            // --- Your coordinate boundary checks & RNTuple filling code here ---
            // Fill your RNTuple branch variables here, including:
            // - Interaction-level flags: var_fiducial, var_mu_length, etc.
            // - Particle-level kinematics: part.end.x, part.pdg, etc.
          
          // 4. WRITE TO RNTUPLE (Writes one entry per interaction containing all its particles)
        
        writer->Fill();
        j++;
      }
    }
    fi++;
    if(fi >= 10 ){break;};
  }
//writer.reset();
  

  std::cout << "SUCCESS: All data successfully written to:" << std::endl;
  std::cout << "output_optimization_ntuple_2.root" << std::endl;

  return 0;

}


