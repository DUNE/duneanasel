///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This header defines the boolean-returning functions to be used in the Nu_muCC0Pi reco interaction selection for ND-LAr

//It contains functions for cut flow LVL2, 3, 4 . LVL1 (Fiducial volume) function can be found in NDLArGeometry.h

//This header also contains functions acting on truth interactions, which are then used in a Var definition
//The Var can be called to identify signal and background sources in the Nu_muCC0Pi selection, via truth-matching


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef NU_MU_CC0PI_FUNCTIONS_H
#define NU_MU_CC0PI_FUNCTIONS_H

//#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
//#include "duneanaobj/StandardRecord"
#include "NDLAr_Geometry.h"
//#include "CAFAna/Core/Var.h"

namespace ana{

namespace NDLArFunctions{

//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//reco selection functions from level 2 to level 4
// you can find level 1 (FV) cut in the geometry header

//reco selection level 2 criteria - intercation produces produces muon
bool ExactlyOnePrimaryMuon(const caf::SRInteractionProxy* sr)
{
    int nMuons = 0;

    for (const auto& p : sr->part.dlp)
    {
        if (p.primary == 1 && p.pdg == 13)
            ++nMuons;
    }

    return nMuons == 1;
}


// reco election level 3 criteria - produced muon has minimum track 20cm
bool PrimMuMinTrackLength(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {
      if (p.primary == 1 && p.pdg == 13 )
      {
        //define track length in here
          double trackstartX = p.start.x;
          double trackstartY = p.start.y;
          double trackstartZ = p.start.z;

          double trackendX = p.end.x;
          double trackendY = p.end.y;
          double trackendZ = p.end.z; 

          double tracklength = std::sqrt(std::pow((trackendX - trackstartX),2) + std::pow((trackendY - trackstartY),2) + std::pow((trackendZ - trackstartZ),2));
        
          if (tracklength >= 20.0){
            return true;
          }
          else{
            return false;
          }
      }
       
    }
    return false;

}

// reco selection level 4 criteria 

bool ProducesZeroPions(const caf::SRInteractionProxy* sr)
{
    int nPions = 0;

    for (const auto& p : sr->part.dlp)
    {
        if (p.primary == 1 && p.pdg == 211)
            ++nPions;
    }

    return nPions == 0;
}


/*
bool RecoParticleContainment(const caf::SRInteractionProxy* sr)
{   //define crop for containment 
    double outer_x_high = 10.0;
    double outer_x_low = 10.0;
    double outer_y_high = 10.0;
    double outer_y_low = 10.0;  
    double outer_z_high = 10.0;
    double outer_z_low = 10.0;

    for (const auto& p : sr->part.dlp){

        double trackstartX = p.start.x;
        double trackstartY = p.start.y;
        double trackstartZ = p.start.z;

        double trackendX = p.end.x;
        double trackendY = p.end.y;
        double trackendZ = p.end.z; 
        
        if (p.primary == 1 && p.pdg == 13){

            //use outerfiducial volume function to check containment in X and Y, but not Z for primary muons
            return NDLArGeo::IsInOuterFiducialVolume(trackstartX, trackstartY, trackstartZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  0.0, outer_z_low) && 
                   NDLArGeo::IsInOuterFiducialVolume(trackendX, trackendY, trackendZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  0.0, outer_z_low);
        }

        else{
            //use outerfiducial volume function to check containment in X, Y and Z for all other particles
            return NDLArGeo::IsInOuterFiducialVolume(trackstartX, trackstartY, trackstartZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  outer_z_high, outer_z_low) && 
                   NDLArGeo::IsInOuterFiducialVolume(trackendX, trackendY, trackendZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  outer_z_high, outer_z_low);
        }
    }
}
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Containemnt:

//function for checking original object type:
std::string PrintRecoObjType(int t){
    //caf::SRRecoBaseID::SRRecoBaseCollectionType t){
  switch(t){
    case 1://caf::SRRecoBaseID::kNDLArDLPShower:
      return "kNDLArDLPShower";
    case 2://caf::SRRecoBaseID::kNDLArDLPTrack:
      return "kNDLArDLPTrack";
    default:
      return std::string("Unknown:") + std::to_string(int(t));
  }
}




// Attempting to create a function for containment:
// loop over each interaction:
//       loop over each particle in the interaction:
//          contain in X
//          contain in Y
//          check if the particle is a primary muon (pdg = 13, primary = 1):
//          if yes, continue (no downstream containment )
//          if no, contain in Z



enum containment_category {
    kcont, kuncontained_back, kuncontained};


containment_category Containment(const caf::SRRecoParticleProxy& p){
    double outer_x_high = 10.0;
    double outer_x_low = 10.0;
    double outer_y_high = 10.0;
    double outer_y_low = 10.0; 
    double outer_z_high = 10.0;
    double outer_z_low = 10.0;

    double trackstartX = p.start.x;
    double trackstartY = p.start.y;
    double trackstartZ = p.start.z;

    double trackendX = p.end.x;
    double trackendY = p.end.y;
    double trackendZ = p.end.z; 

    bool ndlarcont =  NDLArGeo::IsInOuterFiducialVolume(trackstartX, trackstartY, trackstartZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  outer_z_high, outer_z_low) && 
                      NDLArGeo::IsInOuterFiducialVolume(trackendX, trackendY, trackendZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  outer_z_high, outer_z_low);
    if(ndlarcont){
        return kcont;
    } 
    bool cont_with_zero_backface_cut =  NDLArGeo::IsInOuterFiducialVolume(trackstartX, trackstartY, trackstartZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  -100.0, outer_z_low) && 
                   NDLArGeo::IsInOuterFiducialVolume(trackendX, trackendY, trackendZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  -100.0, outer_z_low);
    if(cont_with_zero_backface_cut){
    return kuncontained_back;
    }
    return kuncontained;
    
}


//Original function for containment which at the moment is failing

bool NumuCC0PiPartContainment(const caf::SRInteractionProxy* sr)
{   //define crop for containment 
   
    for (const auto& p : sr->part.dlp){

        auto cont = Containment(p);
        if (p.primary == 1 && p.pdg == 13){
            if(cont == kuncontained){
                return false;
            }
        } else if(cont != kcont){
            return false;
        }
    }
    return true;
}
//split into non-muon containment (bottom bit)
//muon containement (back face OR contained)

//function for all non-muon particle containement 
bool NonMuonsContained(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {
        //skip primary muons
        //if(p.primary == 1 && p.pdg == 13){
         //   continue;
       //}

        //if ((p.pdg == 22 || p.pdg == 11) && (p.end.x <= -1000 && p.end.y<= -1000 && p.end.z<= -1000))
        //{continue;}
        if (p.primary == 1 && p.pdg != 13)
        {

            auto cont = Containment(p);
            if (cont != kcont) {
                return false; //non muon particle escaping the volume fails the cut 
            }
         }
    }
    return true;

}

//non muon TRACK type objects contained 
bool NonMuonsTRACKContained(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {
        //skip primary muons
        //if(p.primary == 1 && p.pdg == 13){
         //   continue;
       //}
        //if ((p.pdg == 22 || p.pdg == 11) && (p.end.x <= -1000 && p.end.y<= -1000 && p.end.z<= -1000))
        //{continue;}

        //function that returns string output of "PrintRecoObjType", defined above this 
        //std::string reco_type_str = PrintRecoObjType(p.type);
            //p.recoobj.type.GetValue());
        // if the particle is a primary, is NOT. muon, and is a Track type
        //if (p.primary == 1 && p.pdg != 13 && (reco_type_str.find("kNDLArDLPTrack") != std::string::npos))
        if (p.primary == 1 && p.pdg != 13 && std::abs(p.pdg) != 11 && std::abs(p.pdg) != 22)
        {

            auto cont = Containment(p);
            if (cont != kcont) {
                return false; //non muon particle escaping the volume fails the cut 
            }
         }
    }
    return true;

}

bool NonMuonsSHOWERContained(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {
        //skip primary muons
        //if(p.primary == 1 && p.pdg == 13){
         //   continue;
       //}
        //if ((p.pdg == 22 || p.pdg == 11) && (p.end.x <= -1000 && p.end.y<= -1000 && p.end.z<= -1000))
        //{continue;}

        //function that returns string output of "PrintRecoObjType", defined above this 
        //std::string reco_type_str = PrintRecoObjType(p.type);
            //p.recoobj.type.GetValue());
        // if the particle is a primary, is NOT. muon, and is a Track type
        //if (p.primary == 1 && p.pdg != 13 && (reco_type_str.find("kNDLArDLPShower") != std::string::npos))
        if (p.primary == 1 && p.pdg != 13 && (std::abs(p.pdg) == 11 || std::abs(p.pdg) == 22))
        {

            auto cont = p.contained;
            if (cont != 1) {
                return false; //non muon particle escaping the volume fails the cut 
            }
         }
    }
    return true;

}





//function for checking primary muons are contained BUT can escape DS
bool PrimaryMuonContained(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {
        if(p.primary ==1 && p.pdg == 13)
        {
            auto cont = Containment(p);

            //If the primary muon is completely uncontained, reject
            if(cont == kuncontained){
                return false;
            }
            //kcont or kuncontained_back are both allowed
        }
    }
    return true;
}



bool NoShowersInteractions(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {//p.primary == 1 && p.recoobj.type != caf::SRRecoBaseID::kNDLArDLPTrack
        //if (p.primary == 1 && p.type !=2 ){
        if (p.primary == 1 && (std::abs(p.pdg) == 11 || std::abs(p.pdg) == 22)) {
            return false;
        }
    }
    return true;
}
/*
//change name to reflect what it does
std::array<caf::SRVector3D,2> FindInteractionBoxCorners(const caf::SRInteractionProxy* sr, double Containment_crop)
{

    if (!sr || sr->part.dlp.empty()) {
        return {caf::SRVector3D{0, 0, 0}, caf::SRVector3D{0, 0, 0}};

        // Initialize mins to max possible value, maxs to lowest possible value
    double x_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::lowest();
    double z_min = std::numeric_limits<double>::max();
    double z_max = std::numeric_limits<double>::lowest();
    
    
    for(const auto& p : sr->part.dlp)
    {
        if(p.primary ==1 && p.pdg == 13){
            continue;
        }
        // Compare against both start and end points for X, Y, and Z
        x_min = std::min({x_min, p.start.x, p.end.x});
        x_max = std::max({x_max, p.start.x, p.end.x});

        y_min = std::min({y_min, p.start.y, p.end.y});
        y_max = std::max({y_max, p.start.y, p.end.y});

        z_min = std::min({z_min, p.start.z, p.end.z});
        z_max = std::max({z_max, p.start.z, p.end.z});
    }
    return {caf::SRVector3D{x_min,y_min,z_min},
             caf::SRVector3D{x_max,y_max,z_max}};
    };
}
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TRUTH functions, used for identifying true signal and background contaminants in the selection
//These functions are implamented in the Var definition below


//check true FV position
bool InFVTrue (const caf::SRTrueInteractionProxy* sr){

  //const std::vector<double> kDefaultBoxMargins = {5.0, 5.0, 5.0}; //in cm
  //const std::vector<double> kDefaultAVMargins = {25.0, 25.0, 25.0}; // in cm

    double outer_x_high = 0.0;
    double outer_x_low = 0.0;
    double outer_y_high = 0.0;
    double outer_y_low = 0.0;  
    double outer_z_high = 0.0;
    double outer_z_low = 0.0;

    //TPC uniform crop
    const double kDefaultTPCMargin = 0.0; // in cm


    double x = sr->vtx.x;
    double y = sr->vtx.y;
    double z = sr->vtx.z;

    

  return NDLArGeo::IsInCombinedNDLArFV(x, y, z, outer_x_high, outer_x_low, 
                                                outer_y_high, outer_y_low, 
                                                outer_z_high, outer_z_low, 
                                                kDefaultTPCMargin);

}

//particle containment 
/*
enum containment_category_truth {
    kcont_truth, kuncontained_back_truth, kuncontained_truth};


containment_category_truth Containment_truth(const caf::SRTrueParticleProxy& p){
    double outer_x_high = 10.0;
    double outer_x_low = 10.0;
    double outer_y_high = 10.0;
    double outer_y_low = 10.0; 
    double outer_z_high = 10.0;
    double outer_z_low = 10.0;

    double trackstartX = p.start_pos.x;
    double trackstartY = p.start_pos.y;
    double trackstartZ = p.start_pos.z;

    double trackendX = p.end_pos.x;
    double trackendY = p.end_pos.y;
    double trackendZ = p.end_pos.z; 

    bool ndlarcont =  NDLArGeo::IsInOuterFiducialVolume(trackstartX, trackstartY, trackstartZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  outer_z_high, outer_z_low) && 
                      NDLArGeo::IsInOuterFiducialVolume(trackendX, trackendY, trackendZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  outer_z_high, outer_z_low);
    if(ndlarcont){
        return kcont_truth;
    } 
    bool cont_with_zero_backface_cut =  NDLArGeo::IsInOuterFiducialVolume(trackstartX, trackstartY, trackstartZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  0.0, outer_z_low) && 
                   NDLArGeo::IsInOuterFiducialVolume(trackendX, trackendY, trackendZ, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  0.0, outer_z_low);
    if(cont_with_zero_backface_cut){
    return kuncontained_back_truth;
    }
    return kuncontained_truth;
    
}


bool NumuCC0PiPartContainmentTruth(const caf::SRTrueInteractionProxy* truth)
{   //define crop for containment 
   
    for (const auto& p : truth->prim){

        auto cont = Containment_truth(p);
        if (p.pdg == 13){
            if(cont == kuncontained_truth){
                return false;
            }
        } else if(cont != kcont_truth){
            return false;
        }
    }
    return true;
}



*/

//check true muon track length
bool PrimaryTrueMuonHasMinimumTrackLength(const caf::SRTrueInteractionProxy* truth)
{
    for (const auto& p : truth->prim)
    {
        if (p.pdg == 13)
        {
            const double dx = p.end_pos.x - p.start_pos.x;
            const double dy = p.end_pos.y - p.start_pos.y;
            const double dz = p.end_pos.z - p.start_pos.z;

            const double tracklength = std::sqrt(dx*dx + dy*dy + dz*dz);

            return tracklength >= 20.0;
        }
    }

    return false;
}


//Returns truth matched interaction for a reco interaction input 
const caf::SRTrueInteractionProxy* GetTruthMatch(const caf::SRInteractionProxy* sr)
{
    if (sr->truth.empty())
        return nullptr;

    // Find best-matched true interaction
    size_t tidx = 0;
    float maxOverlap = 0.0;

    for (size_t i = 0; i < sr->truthOverlap.size(); ++i) {
        if (sr->truthOverlap[i] > maxOverlap) {
            maxOverlap = sr->truthOverlap[i];
            tidx = i;
        }
    }

    // Get ancestor SRProxy
    const caf::SRProxy* proxy = sr->Ancestor<caf::SRProxy>();

    // Return matched true interaction
    return &proxy->mc.nu[sr->truth[tidx]];
}

//enumerates the interaction type 
enum TruthInteractionType {
    kNuMuCC0PiFV = 0,
    kNuMuCCnPiFV = 1,
    kNuMu_NuMuBarCCincOOFV = 2,
    kNCinc       = 3,
    kNuMuBarCCincFV = 4,
    kNuECCinc    = 5,
    kNuEBarCCinc = 6,
    kNoMatch     = 7,
    kUncatagorisable = 8,
    kNTypes = 9
};

std::string TruthInteractionTypeToString(TruthInteractionType tit){
    switch(tit){
        case kNuMuCC0PiFV: return "NuMuCC0PiFV";
        case kNuMuCCnPiFV: return "NuMuCCnPiFV";
        case kNuMu_NuMuBarCCincOOFV: return "NuMu_NuMuBarCCincOOFV";
        case kNCinc: return "NCinc";
        case kNuMuBarCCincFV: return "NuMuBarCCincFV";
        case kNuECCinc: return "NuECCinc";
        case kNuEBarCCinc: return "NuEBarCCinc";
        case kNoMatch: return "NoMatch";
        case kUncatagorisable: return "Uncatagorisable";
        default: return "Unknown";
    }
}
std::string TruthInteractionTypeToString(int tit){
    return TruthInteractionTypeToString(TruthInteractionType(tit));
}








//generate Var enumerates to correct interaction type 
inline double GetTruthInteractionType(const caf::SRInteractionProxy* sr){ 
    const caf::SRTrueInteractionProxy* truthInt = GetTruthMatch(sr);
    if (!truthInt) 
        return static_cast<double>(kNoMatch);

    // ------------------
    //basic truth info
    //-------------------

    const int nupdg = truthInt->pdg;
    const bool nu_iscc = truthInt->iscc;

    const bool inFV =
            InFVTrue(truthInt);
    
  
    const bool zeroPi = (truthInt->npip == 0 && truthInt->npim == 0 && truthInt->npi0 == 0);

    const bool hasLongMuon = PrimaryTrueMuonHasMinimumTrackLength(truthInt);

    const bool nPi = (truthInt->npip > 0 || truthInt->npim > 0 || truthInt->npi0 > 0);


    //------------------------------
    //NC check.
    //------------------------------

    if (!nu_iscc) {
        return static_cast<double>(kNCinc);
    }

    //------------------------------
    // nu_mu CC check
    //------------------------------

    else if (nupdg == 14) {


            //nu_mu CC interaction outside the FV
            if(!inFV && hasLongMuon){
                return static_cast<double>(kNuMu_NuMuBarCCincOOFV);
            }
            // nu_mu CC 0Pi with a truth muon track >= 20 cm AND particles contained/ muons DS
            if(zeroPi && hasLongMuon){
                return static_cast<double>(kNuMuCC0PiFV);
                
            
            }
            //nu_mu CC with at least one pions
            else if(nPi){

                return static_cast<double>( kNuMuCCnPiFV);
            }

            else {
                //nu_mu CC in FV by not nPi
                //this catches 0Pi as well as anything unusual
                return static_cast<double>(kUncatagorisable);
            }

    }

    //------------------------------
    //anti-nu_mu cc
    //------------------------------

    else if (nupdg == -14) {

            if (!inFV && hasLongMuon) {
                return static_cast<double>(kNuMu_NuMuBarCCincOOFV);
            }
        
        //anti-nu_mu CC inside FV
        return static_cast<double >(kNuMuBarCCincFV);


    }

    //------------------------------
    //nu_e CC check
    //------------------------------


    else if (nupdg == 12) {
            return static_cast<double>(kNuECCinc);
    }   


    //------------------------------
    //anti-nu_e CC check
    //------------------------------

    else if(nupdg == -12){
            return static_cast<double>(kNuEBarCCinc);
    }

    //-----------------------------
    //if we get here, we have a CC interaction that is not nu_mu or nu
    //-----------------------------
    else{
        return static_cast<double>(kUncatagorisable);
    }

}


    




} // namespace NDLArFunctions
} // namespace ana
#endif // NU_MU_CC0PI_FUNCTIONS_H



