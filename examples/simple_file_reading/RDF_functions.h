#include <ROOT/RVec.hxx>
#include "duneanaobj/StandardRecord/SRTrueInteraction.h"
#include "duneanaobj/StandardRecord/SRTrueParticle.h"
#include "duneanaobj/StandardRecord/SRVector3D.h"
#include "duneanaobj/StandardRecord/SRLorentzVector.h"
ROOT::VecOps::RVec<Float_t> nu_prop_by_pdg(const ROOT::VecOps::RVec<Float_t>& p, const ROOT::VecOps::RVec<Int_t>& pdg, Int_t pdg_target) {

    ROOT::VecOps::RVec<Float_t> OUT;

    for (size_t i = 0; i < p.size(); ++i) {
        if (pdg[i] == pdg_target) {
            OUT.push_back(p[i]);
        }
    }

    return OUT;
}

ROOT::VecOps::RVec<Int_t> nu_prop_by_pdg(const ROOT::VecOps::RVec<Int_t>& p, const ROOT::VecOps::RVec<Int_t>& pdg, Int_t pdg_target) {

    ROOT::VecOps::RVec<Int_t> OUT;

    for (size_t i = 0; i < p.size(); ++i) {
        if (pdg[i] == pdg_target) {
            OUT.push_back(p[i]);
        }
    }

    return OUT;
}

ROOT::VecOps::RVec<caf::SRLorentzVector> nu_prop_by_pdg(const ROOT::VecOps::RVec<caf::SRLorentzVector>& p, const ROOT::VecOps::RVec<Int_t>& pdg, Int_t pdg_target) {

    ROOT::VecOps::RVec<caf::SRLorentzVector> OUT;

    for (size_t i = 0; i < p.size(); ++i) {
        if (pdg[i] == pdg_target) {
            OUT.push_back(p[i]);
        }
    }

    return OUT;
}

ROOT::VecOps::RVec<caf::SRVector3D> nu_prop_by_pdg(const ROOT::VecOps::RVec<caf::SRVector3D>& p, const ROOT::VecOps::RVec<Int_t>& pdg, Int_t pdg_target) {

    ROOT::VecOps::RVec<caf::SRVector3D> OUT;

    for (size_t i = 0; i < p.size(); ++i) {
        if (pdg[i] == pdg_target) {
            OUT.push_back(p[i]);
        }
    }

    return OUT;
}

ROOT::VecOps::RVec<Int_t> prim_pdg(const ROOT::VecOps::RVec<std::vector<caf::SRTrueParticle>>& p) {

    ROOT::VecOps::RVec<Int_t> OUT;

    for (const auto & prim_vec : p) {
            for (const auto & prim : prim_vec) {
                OUT.push_back(prim.pdg);
        } 
    }

    return OUT;
}