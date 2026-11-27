#pragma once

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana::numucc0pi {

//enumerates the interaction type
enum TruthInteractionType {
    kNuMuCC0PiFV_ArgTarg = 0,
    kNuMuCCnPiFV = 1,
    kNuMu_NuMuBarCCincOOFV = 2,
    kNCinc       = 3,
    kNuMuBarCCincFV = 4,
    kNuECCinc    = 5,
    kNuEBarCCinc = 6,
    kNoMatch     = 7,
    kUncatagorisable = 8,
    kNuMuCC0PiFV_NONArgTarg = 9,
    kNTypes = 10

};

std::string TruthInteractionTypeToString(TruthInteractionType tit){
    switch(tit){
        case kNuMuCC0PiFV_ArgTarg: return "NuMuCC0Pi_Ar_FV";
        case kNuMuCCnPiFV: return "NuMuCCnPiFV";
        case kNuMu_NuMuBarCCincOOFV: return "NuMu_NuMuBarCCincOOFV";
        case kNCinc: return "NCinc";
        case kNuMuBarCCincFV: return "NuMuBarCCincFV";
        case kNuECCinc: return "NuECCinc";
        case kNuEBarCCinc: return "NuEBarCCinc";
        case kNoMatch: return "NoMatch";
        case kUncatagorisable: return "Uncatagorisable";
        case kNuMuCC0PiFV_NONArgTarg: return "kNuMuCC0Pi_NONAr_FV";
        default: return "Unknown";
    }
}
std::string TruthInteractionTypeToString(int tit){
    return TruthInteractionTypeToString(TruthInteractionType(tit));
}

//Returns a true interaction's interaction type
inline double TruthInteractionChecker(const caf::SRTrueInteractionProxy* truthInt){
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

        //const bool hasLongMuon = PrimaryTrueMuonHasMinimumTrackLength(truthInt);

        const bool nPi = (truthInt->npip > 0 || truthInt->npim > 0 || truthInt->npi0 > 0);

        const bool Is_Arg = truthInt->targetPDG == 1000180400;

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

                //ALTERNATIVE FOR Nu_muCCINC selection, no pion selection happening here
                /*
                if(!inFV){
                    return static_cast<double>(kNuMu_NuMuBarCCincOOFV);
                }
                else{
                // nu_mu CC 0Pi with a truth muon track >= 20 cm AND particles contained/ muons DS
            //else if(zeroPi ){
                    // Split based on whether the struck target is Argon-40 or not
                    if (Is_Arg) {
                        return static_cast<double>(kNuMuCC0PiFV_ArgTarg);
                    } else {
                        return static_cast<double>(kNuMuCC0PiFV_NONArgTarg);
                    }
                }
                */

                //nu_mu CC interaction outside the FV

            if(!inFV ){
                return static_cast<double>(kNuMu_NuMuBarCCincOOFV);
            }

            // nu_mu CC 0Pi with a truth muon track >= 20 cm AND particles contained/ muons DS
            else if(zeroPi ){
            // Split based on whether the struck target is Argon-40 or not
                if (Is_Arg) {
                    return static_cast<double>(kNuMuCC0PiFV_ArgTarg);
                } else {
                    return static_cast<double>(kNuMuCC0PiFV_NONArgTarg);
                }
            }
            //nu_mu CC with at least one pions
            else if(nPi){
                return static_cast<double>( kNuMuCCnPiFV);
            }
            else {
                //nu_mu CC in FV by not nPi
                // //this catches 0Pi as well as anything unusual
                return static_cast<double>(kUncatagorisable);}
        }


            /*
                //nu_mu CC with at least one pions
                else if(nPi){

                    return static_cast<double>( kNuMuCCnPiFV);
                }

                else {
                    //nu_mu CC in FV by not nPi
                    //this catches 0Pi as well as anything unusual
                    return static_cast<double>(kUncatagorisable);
                }
    */


        //------------------------------
        //anti-nu_mu cc
        //------------------------------

        else if (nupdg == -14) {

                if (!inFV /*&& hasLongMuon*/) {
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

    //Takes a reco intercation, finds its best-matched truth interaction, and returns the interaction type.
    inline double GetTruthInteractionType(const caf::SRInteractionProxy* sr){
        const caf::SRTrueInteractionProxy* truthInt = GetTruthMatch(sr);

        return TruthInteractionChecker(truthInt);

    }


}

 // namespace ana::numucc0pi
