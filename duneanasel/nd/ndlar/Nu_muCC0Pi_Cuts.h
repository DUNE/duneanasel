///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This header contains Cut definitions to be used in the Nu_muCC0Pi reco selection cut flow

//The fiducial volume (LVL1) cut uses the FV function from NDLArGeometry.h

//the other cuts use boolean-returning functions defined in Nu_muCC0Pi_functions.h

//These cuts can be implamented in the CAFAna framework for the Nu_muCC0Pi selection

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ND_LAR_CUTS_H
#define ND_LAR_CUTS_H

#include "CAFAna/Core/Cut.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanasel/nd/ndlar/NDLAr_Geometry.h"
#include "duneanasel/nd/ndlar/Nu_muCC0Pi.h"
#include <vector>

namespace ana {

//selection level 1 criteria - purely FV cut only 
const Cut kNDLArFiducialVolume([](const caf::SRInteractionProxy* sr) {
    //const std::vector<double> kDefaultBoxMargins = {5.0, 5.0, 5.0};  // in cm
    //const std::vector<double> kDefaultAVMargins  = {50.0, 50.0, 50.0}; // in cm

    //Outer FV crops
    double outer_x_high = 50.0;
    double outer_x_low = 50.0;
    double outer_y_high = 50.0;
    double outer_y_low = 50.0;  
    double outer_z_high = 150.0;
    double outer_z_low = 50.0;
    
    //TPC uniform crop
    const double kDefaultTPCMargin = 5.0; // in cm

    double x = sr->vtx.x;
    double y = sr->vtx.y;
    double z = sr->vtx.z;

    return NDLArGeo::IsInCombinedNDLArFV(x, y, z, outer_x_high, outer_x_low, 
                                                  outer_y_high, outer_y_low, 
                                                  outer_z_high, outer_z_low, 
                                                  kDefaultTPCMargin );
});

//Selection level 1.5 criteria - particle containment 
const Cut kNDLArPartContainment (NDLArFunctions::NumuCC0PiPartContainment);

//selection level 2 criteria - require 1 primary muon produced by interaction
const Cut kNDLArPrimMu (NDLArFunctions::ExactlyOnePrimaryMuon);

//selection level 3 criteria - require primary muon track length to be 20cm minimum
const Cut kNDLArMuTrackLength(NDLArFunctions::PrimMuMinTrackLength);

//selection level 4 criteria - require no pions produced (soon to be renamed lvl 5)
const Cut kNDLArZeroPions(NDLArFunctions::ProducesZeroPions);

//selection level 5 criteria - require all particles to be contained
//const Cut kNDLArAllParticlesContained(NDLArFunctions::NumuCC0PiPartContainment);]

//new selection level 5 criteria - require all primary muons to be contained, but can escape DS
const Cut kMuonsContained(NDLArFunctions::PrimaryMuonContained);

//new selection level 6 criteria - require ALL non-muon TRACK particles to be contained
const Cut kNonMuonsTRACKContained(NDLArFunctions::NonMuonsTRACKContained);

//new selection level 7 criteria - require ALL non-muon SHOWERS particles to be contained
const Cut kNonMuonsSHOWERContained(NDLArFunctions::NonMuonsSHOWERContained);

//new selection criteria - if an interaction produces a shower, reject it 
const Cut kNonShowerInteractions(NDLArFunctions::NoShowersInteractions);

/*
//new containment cut which takes crop as an input parameter
const Cut NonMuonTracksCont_BoxMethod(NDLArFunctions::NonMuonTracksCont_BoxMethod);
*/

} // namespace ana

#endif // ND_LAR_CUTS_H