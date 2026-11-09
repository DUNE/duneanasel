#include "duneanaobj/StandardRecord/SRInteraction.h"

namespace ana {

constexpr static auto kCannotFindTruth = std::numeric_limits<size_t>::max();

inline size_t GetTrueInteractionIndex(caf::SRInteraction const &nu_int) {
  if (!nu_int.truthOverlap.size()) {
    return kCannotFindTruth;
  }
  size_t idx = std::distance(
      nu_int.truthOverlap.begin(),
      std::max_element(nu_int.truthOverlap.begin(), nu_int.truthOverlap.end()));
  if (idx >= nu_int.truth.size()) {
    return kCannotFindTruth;
  }
  return nu_int.truth[idx];
}



//////////////////////////////////////////////////////////////////////////////////////////////




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
        default: return "Unknownget";
    }
}
std::string TruthInteractionTypeToString(int tit){
    return TruthInteractionTypeToString(TruthInteractionType(tit));
}









inline int GetTruthInteractionType(const caf::SRInteractionProxy* sr){ 
    const caf::SRTrueInteractionProxy* truthInt = GetTruthMatch(sr);
    if (!truthInt) 
        return static_cast<double>(kNoMatch);

    // ------------------
    //basic truth info
    //-------------------

    const int nupdg = truthInt->pdg;
    const bool nu_iscc = truthInt->iscc;

    const bool inFV =
            ana::NDLArFunctions::InFVTrue(truthInt);
    
  
    const bool zeroPi = (truthInt->npip == 0 && truthInt->npim == 0 && truthInt->npi0 == 0);

    const bool hasLongMuon = ana::NDLArFunctions::PrimaryTrueMuonHasMinimumTrackLength(truthInt);

    const bool nPi = (truthInt->npip > 0 || truthInt->npim > 0 || truthInt->npi0 > 0);


    //------------------------------
    //NC check.
    //------------------------------

    if (!nu_iscc) {
        return kNCinc;
    }

    //------------------------------
    // nu_mu CC check
    //------------------------------

    else if (nupdg == 14) {


            //nu_mu CC interaction outside the FV
            if(!inFV && hasLongMuon){
                return kNuMu_NuMuBarCCincOOFV;
            }
            // nu_mu CC 0Pi with a truth muon track >= 20 cm AND particles contained/ muons DS
            if(zeroPi && hasLongMuon){
                return kNuMuCC0PiFV;
                
            
            }
            //nu_mu CC with at least one pions
            else if(nPi){

                return  kNuMuCCnPiFV;
            }

            else {
                //nu_mu CC in FV by not nPi
                //this catches 0Pi as well as anything unusual
                return kUncatagorisable;
            }

    }

    //------------------------------
    //anti-nu_mu cc
    //------------------------------

    else if (nupdg == -14) {

            if (!inFV && hasLongMuon) {
                return kNuMu_NuMuBarCCincOOFV;
            }
        
        //anti-nu_mu CC inside FV
        return kNuMuBarCCincFV;


    }

    //------------------------------
    //nu_e CC check
    //------------------------------


    else if (nupdg == 12) {
            return kNuECCinc;
    }   


    //------------------------------
    //anti-nu_e CC check
    //------------------------------

    else if(nupdg == -12){
            return kNuEBarCCinc;
    }

    //-----------------------------
    //if we get here, we have a CC interaction that is not nu_mu or nu
    //-----------------------------
    else{
        return kUncatagorisable;
    }

}


    





} // namespace ana