///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This header defines the boolean-returning functions to be used in the Nu_muCC0Pi reco interaction selection for ND-LAr
//For selection functions with variable parameters, a version of the function has been incldued which is optimised for a Nu_muCC0Pi selection
// This header  also contains helper functions  


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "duneanasel/nd/ndlar/NDLAr_Geometry.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

#include "duneanaobj/StandardRecord/SREnums.h"

#include "Math/Point3D.h"
#include <limits>
#include <algorithm>
#include <cmath>

namespace ana::numucc0pi {


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions used for selection in Nu_mu CC0Pi selection within duneanasel:

//Primary-shower-producing interaction removal:

// NO interactions with primary Showers.                       
bool NoShowersInteractions(const caf::SRInteractionProxy* sr) 
{
    for (const auto& p : sr->part.dlp)
    {
        if (p.recoobj.type == caf::SRRecoBaseID::SRRecoBaseCollectionType:: kNDLArDLPShower){
            return false;
        }
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//fiducial volume, with input parameters (using detector geometry from NDLAr_Geometry )

bool IsInNDLArFiducialVolume_variable(const caf::SRInteractionProxy* sr,
                                      double X_crop,
                                      double Y_crop,
                                      double Z_low_crop,
                                      double Z_high_crop,
                                      double TPC_crop)
{
    double reco_vtx_x = sr->vtx.x;
    double reco_vtx_y = sr->vtx.y;
    double reco_vtx_z = sr->vtx.z;

    return ana::NDLArGeo::IsInCombinedNDLArFV(reco_vtx_x, reco_vtx_y, reco_vtx_z,
                                              X_crop, X_crop,
                                              Y_crop, Y_crop,
                                              Z_high_crop, Z_low_crop,
                                              TPC_crop);
}

//level 2 - OPTIMIZED fiducial volume  (using detector geometry from NDLAr_Geometry )

bool IsInNDLArFiducialVolume(const caf::SRInteractionProxy* sr) 
{

double X_crop =10 ;
double Y_crop = 10 ;
double Z_low_crop = 40 ; 
double Z_high_crop = 20 ; 
double TPC_crop = 2.0 ; 


return IsInNDLArFiducialVolume_variable(sr,
                                      X_crop,
                                      Y_crop,
                                      Z_low_crop,
                                      Z_high_crop,
                                      TPC_crop);

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Nu_muCC selection, intercation produces produces muon 
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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Track length checker - Helper function

double PartTrackLength(const caf::SRRecoParticleProxy& p)
{
    double trackstartX = p.start.x;
    double trackstartY = p.start.y;
    double trackstartZ = p.start.z;
    double trackendX = p.end.x;
    double trackendY = p.end.y;
    double trackendZ = p.end.z;

    double tracklength = std::sqrt(std::pow((trackendX - trackstartX),2) + std::pow((trackendY - trackstartY),2) + std::pow((trackendZ - trackstartZ),2));
    return tracklength;
}

// muon has minimum track length - with variable parameter for track_minimum 
bool PrimMuMinTrackLength_variable(const caf::SRInteractionProxy* sr, 
                                    double track_min)
{
    for (const auto& p : sr->part.dlp)
    {
      if (p.primary == 1 && p.pdg == 13 )
      {
          if (PartTrackLength(p) >= track_min){
            return true;
          }
          else{
            return false;
          }
      }
    }
    return false;
}


// produced muon has OPTIMISED minimum track 20cm
bool PrimMuMinTrackLength(const caf::SRInteractionProxy* sr)
{
    return PrimMuMinTrackLength_variable(sr, 18.0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Interaction produces zero pions of neutral, +, - charge. 

bool ProducesZeroPions(const caf::SRInteractionProxy* sr) 
{
    int nPions = 0;
    for (const auto& p : sr->part.dlp)
    {
        if (p.primary == 1 && std::abs(p.pdg) == 211)
            ++nPions;
    }
    return nPions == 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//helper function for checking that a box is within a box, used for containment 

inline bool BoxContainedInBox( double x_min, double y_min, double z_min,
                                  double x_max, double y_max, double z_max,  double Containment_crop)
    {
        return (ana::NDLArGeo::IsInOuterFiducialVolume(x_min,y_min,z_min,Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop ) 
            &&ana::NDLArGeo::IsInOuterFiducialVolume(x_max,y_max,z_max,  Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop ));
    }




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//level 6 Containment VARIABLE split contianment functions: non-muon track containment and muon track containment 
//containment using a bounding box from muon particle tracks 

// Checks containment of the bounding box built strictly from non-muon tracks (all 6 faces)
bool NonMuonTracksContained_variable(const caf::SRInteractionProxy* sr, double Containment_Crop = 5.0)
{
    if (!sr) return false;

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    double zmin = std::numeric_limits<double>::max();
    double zmax = std::numeric_limits<double>::lowest();
    int non_muon_track_count = 0;

    for (const auto& part : sr->part.dlp) {
        if (part.primary != 1) continue;
        bool is_primary_muon = (part.primary == 1 && std::abs(part.pdg) == 13);
        bool is_shower = (part.recoobj.type == caf::SRRecoBaseID::SRRecoBaseCollectionType::kNDLArDLPShower);
        bool is_track = !is_shower;

        if (is_track && !is_primary_muon) {
            non_muon_track_count++;

            xmin = std::min({xmin, static_cast<double>(part.start.x), static_cast<double>(part.end.x)});
            xmax = std::max({xmax, static_cast<double>(part.start.x), static_cast<double>(part.end.x)});

            ymin = std::min({ymin, static_cast<double>(part.start.y), static_cast<double>(part.end.y)});
            ymax = std::max({ymax, static_cast<double>(part.start.y), static_cast<double>(part.end.y)});

            zmin = std::min({zmin, static_cast<double>(part.start.z), static_cast<double>(part.end.z)});
            zmax = std::max({zmax, static_cast<double>(part.start.z), static_cast<double>(part.end.z)});
        }
    }

    // Fallback to vertex position if no non-muon tracks exist
    if (non_muon_track_count == 0) {
        xmin = xmax = sr->vtx.x;
        ymin = ymax = sr->vtx.y;
        zmin = zmax = sr->vtx.z;
    }

    return ana::NDLArGeo::BoxContainedInBox(xmin, ymin, zmin,
                                          xmax, ymax, zmax,
                                          Containment_Crop);
}



bool MuonTracksContained_variable(const caf::SRInteractionProxy* sr, double Containment_Crop = 5.0)
{
    if (!sr) return false;

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::lowest();
    double zmin = std::numeric_limits<double>::max();
    double zmax = std::numeric_limits<double>::lowest();
    int muon_track_count = 0;

    for (const auto& part : sr->part.dlp) {
        bool is_primary_muon = (part.primary == 1 && std::abs(part.pdg) == 13);
        bool is_shower = (part.recoobj.type == caf::SRRecoBaseID::SRRecoBaseCollectionType::kNDLArDLPShower);
        bool is_track = !is_shower;

        if (is_track && is_primary_muon) {
            muon_track_count++;

            xmin = std::min({xmin, static_cast<double>(part.start.x), static_cast<double>(part.end.x)});
            xmax = std::max({xmax, static_cast<double>(part.start.x), static_cast<double>(part.end.x)});

            ymin = std::min({ymin, static_cast<double>(part.start.y), static_cast<double>(part.end.y)});
            ymax = std::max({ymax, static_cast<double>(part.start.y), static_cast<double>(part.end.y)});

            zmin = std::min({zmin, static_cast<double>(part.start.z), static_cast<double>(part.end.z)});
            zmax = std::max({zmax, static_cast<double>(part.start.z), static_cast<double>(part.end.z)});
        }
    }

    // Fallback to vertex position if no muon tracks exist
    if (muon_track_count == 0) {
        xmin = xmax = sr->vtx.x;
        ymin = ymax = sr->vtx.y;
        zmin = zmax = sr->vtx.z;
    }

    const double x_lo = ana::NDLArGeo::NDLArXLo + Containment_Crop;
    const double x_hi = ana::NDLArGeo::NDLArXHi - Containment_Crop;
    const double y_lo = ana::NDLArGeo::NDLArYLo + Containment_Crop;
    const double y_hi = ana::NDLArGeo::NDLArYHi - Containment_Crop;
    const double z_lo = ana::NDLArGeo::NDLArZLo + Containment_Crop;

    return (xmin >= x_lo && xmax <= x_hi &&
            ymin >= y_lo && ymax <= y_hi &&
            zmin >= z_lo);
}
    
bool Nu_muCCContained(const caf::SRInteractionProxy* sr,double Containment_Crop = 5.0)
{
    return NonMuonTracksContained_variable(sr, Containment_Crop) && MuonTracksContained_variable(sr, Containment_Crop); 
    
    
}


} // namespace ana
