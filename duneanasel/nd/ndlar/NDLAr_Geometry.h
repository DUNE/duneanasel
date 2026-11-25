
#ifndef NDLARGEOMETRY_H
#define NDLARGEOMETRY_H

#include <vector>
#include <cmath>
#include "Math/Point3D.h"
#include <limits>
#include <algorithm>


namespace ana{

namespace NDLArGeo {

    // Structure defining the bounding box of an individual TPC volume
    struct Box {
        double x_min, x_max;
        double y_min, y_max;
        double z_min, z_max;

        bool Contains(double x, double y, double z) const {
            return (x >= x_min && x <= x_max &&
                    y >= y_min && y <= y_max &&
                    z >= z_min && z <= z_max);
        }
    };


    // define the full TPC Dimensions from YAML file
    const double TPC_DX = 46.788;
    const double TPC_DY = 297.6;
    const double TPC_DZ = 95.232;

    const double HALF_DX = TPC_DX / 2.0;
    const double HALF_DY = TPC_DY / 2.0;
    const double HALF_DZ = TPC_DZ / 2.0;


 //new definitions for detector geometry, 
 //using maximum boundaries of TPC modules
    const double NDLArXLo = -323.7115 - HALF_DX;
    const double NDLArXHi = 323.7115 + HALF_DX;
    const double NDLArYLo = -66.8713 - HALF_DY;
    const double NDLArYHi = -66.8713 + HALF_DY;
    const double NDLArZLo = 465.7559 - HALF_DZ;
    const double NDLArZHi = 865.7559 + HALF_DZ;

    // Center positions of all 70 TPCs [X, Y, Z]
    inline const std::vector<std::vector<double>> TPC_CENTERS = {
        {-276.2885, -66.8713, 465.7559}, {-323.7115, -66.8713, 465.7559},
        {-276.2885, -66.8713, 565.7559}, {-323.7115, -66.8713, 565.7559},
        {-276.2885, -66.8713, 665.7559}, {-323.7115, -66.8713, 665.7559},
        {-276.2885, -66.8713, 765.7559}, {-323.7115, -66.8713, 765.7559},
        {-276.2885, -66.8713, 865.7559}, {-323.7115, -66.8713, 865.7559},
        {-176.2885, -66.8713, 465.7559}, {-223.7115, -66.8713, 465.7559},
        {-176.2885, -66.8713, 565.7559}, {-223.7115, -66.8713, 565.7559},
        {-176.2885, -66.8713, 665.7559}, {-223.7115, -66.8713, 665.7559},
        {-176.2885, -66.8713, 765.7559}, {-223.7115, -66.8713, 765.7559},
        {-176.2885, -66.8713, 865.7559}, {-223.7115, -66.8713, 865.7559},
        { -76.2885, -66.8713, 465.7559}, {-123.7115, -66.8713, 465.7559},
        { -76.2885, -66.8713, 565.7559}, {-123.7115, -66.8713, 565.7559},
        { -76.2885, -66.8713, 665.7559}, {-123.7115, -66.8713, 665.7559},
        { -76.2885, -66.8713, 765.7559}, {-123.7115, -66.8713, 765.7559},
        { -76.2885, -66.8713, 865.7559}, {-123.7115, -66.8713, 865.7559},
        {  23.7115, -66.8713, 465.7559}, { -23.7115, -66.8713, 465.7559},
        {  23.7115, -66.8713, 565.7559}, { -23.7115, -66.8713, 565.7559},
        {  23.7115, -66.8713, 665.7559}, { -23.7115, -66.8713, 665.7559},
        {  23.7115, -66.8713, 765.7559}, { -23.7115, -66.8713, 765.7559},
        {  23.7115, -66.8713, 865.7559}, { -23.7115, -66.8713, 865.7559},
        { 123.7115, -66.8713, 465.7559}, {  76.2885, -66.8713, 465.7559},
        { 123.7115, -66.8713, 565.7559}, {  76.2885, -66.8713, 565.7559},
        { 123.7115, -66.8713, 665.7559}, {  76.2885, -66.8713, 665.7559},
        { 123.7115, -66.8713, 765.7559}, {  76.2885, -66.8713, 765.7559},
        { 123.7115, -66.8713, 865.7559}, {  76.2885, -66.8713, 865.7559},
        { 223.7115, -66.8713, 465.7559}, { 176.2885, -66.8713, 465.7559},
        { 223.7115, -66.8713, 565.7559}, { 176.2885, -66.8713, 565.7559},
        { 223.7115, -66.8713, 665.7559}, { 176.2885, -66.8713, 665.7559},
        { 223.7115, -66.8713, 765.7559}, { 176.2885, -66.8713, 765.7559},
        { 223.7115, -66.8713, 865.7559}, { 176.2885, -66.8713, 865.7559},
        { 323.7115, -66.8713, 465.7559}, { 276.2885, -66.8713, 465.7559},
        { 323.7115, -66.8713, 565.7559}, { 276.2885, -66.8713, 565.7559},
        { 323.7115, -66.8713, 665.7559}, { 276.2885, -66.8713, 665.7559},
        { 323.7115, -66.8713, 765.7559}, { 276.2885, -66.8713, 765.7559},
        { 323.7115, -66.8713, 865.7559}, { 276.2885, -66.8713, 865.7559}
    };

    // Structure describing an axis-aligned bounding box
    struct AABBCrop {
        double xmin, xmax, ymin, ymax, zmin, zmax;
    };

    AABBCrop MakeUniformCrop(double const &c) { return AABBCrop{c, c, c, c, c, c}; }

    struct AABB {
        ROOT::Math::XYZPoint min, max;

        bool CroppedContains(double const &x, double const &y, double const &z,
                            AABBCrop const &crop) const {
            return ((x >= (min.x() + crop.xmin)) && (x <= (max.x() - crop.xmax)) &&
                    (y >= (min.y() + crop.ymin)) && (y <= (max.y() - crop.ymax)) &&
                    (z >= (min.z() + crop.zmin)) && (z <= (max.z() - crop.zmax)));
        }

        bool CroppedContains(ROOT::Math::XYZPoint const &p,
                            AABBCrop const &crop) const {
            return CroppedContains(p.x(), p.y(), p.z(), crop);
        }

        bool Contains(double const &x, double const &y, double const &z) const {
            return ((x >= min.x()) && (x <= max.x()) && (y >= min.y()) &&
                    (y <= max.y()) && (z >= min.z()) && (z <= max.z()));
        }

        bool Contains(ROOT::Math::XYZPoint const &p) const {
            return Contains(p.x(), p.y(), p.z());
        }
    };

    inline auto GetTPCAABBs() {
        std::vector<AABB> tpcs;

        for (const auto &center : TPC_CENTERS) {
            tpcs.push_back(
                AABB{{center[0] - HALF_DX, center[1] - HALF_DY, center[2] - HALF_DZ},
                    {center[0] + HALF_DX, center[1] + HALF_DY, center[2] + HALF_DZ}});
        }
        return tpcs;
    }

    inline auto GetNDLArAABB() {
        return AABB{{NDLArXLo, NDLArYLo, NDLArZLo}, {NDLArXHi, NDLArYHi, NDLArZHi}};
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Position Checking functions:

    // Checks if position (x, y, z) is inside the outer active volume box
    inline bool IsInOuterFiducialVolume(double x, double y, double z, 
                                        double margin_x_high = 0.0, double margin_x_low = 0.0,
                                        double margin_y_high=0.0, double margin_y_low=0.0,
                                        double margin_z_high = 0.0, double margin_z_low = 0.0){
                                        
        static auto const NDLArAABB = GetNDLArAABB();

        return NDLArAABB.CroppedContains(x, y, z,
                                   AABBCrop{margin_x_low, margin_x_high,
                                            margin_y_low, margin_y_high,
                                            margin_z_low, margin_z_high});
    }

    // Checks if position (x, y, z) is inside any of the individual TPC boxes
    inline bool IsInBoxesFiducialVolume(double x, double y, double z,
                                        double margin_x = 0.0,
                                        double margin_y = 0.0,
                                        double margin_z = 0.0) {
        static auto const TPCAABBs = GetTPCAABBs();

        for (const auto &tpc : TPCAABBs) {
            if (tpc.CroppedContains(x, y, z,
                                    AABBCrop{margin_x, margin_x, margin_y, margin_y,
                                            margin_z, margin_z})) {
            return true;
            }
        }
        return false;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Function which takes the following inputs:
    //vertex position x, y, z
    //TPC box FV margins x, y, z
    //Outer FV margins x, y, z, both high and low

    inline bool IsInCombinedNDLArFV(double x, double y, double z, 
                                    double o_mx_high, double o_mx_low, 
                                    double o_my_high, double o_my_low,
                                    double o_mz_high, double o_mz_low,
                                    double box_margin){
       

        if (!IsInOuterFiducialVolume(x, y, z,o_mx_high ,o_mx_low, o_my_high, o_my_low, o_mz_high, o_mz_low)) return false;

        
        return IsInBoxesFiducialVolume(x, y, z, box_margin, box_margin, box_margin);

    }


    inline bool BoxContainedInBox( double x_min, double y_min, double z_min,
                                  double x_max, double y_max, double z_max,  double Containment_crop)
    {
        return (IsInOuterFiducialVolume(x_min,y_min,z_min,Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop ) 
            &&IsInOuterFiducialVolume(x_max,y_max,z_max,  Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop,
                                                          Containment_crop,Containment_crop ));
    }
    
} //Namespace NDLArGeo
} //Namespace ana
#endif // NDLARGEOMETRY_H
