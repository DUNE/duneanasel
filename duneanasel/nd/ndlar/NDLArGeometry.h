#ifndef NDLARGEOMETRY_H
#define NDLARGEOMETRY_H

#include "Math/Point3D.h"

#include <cmath>
#include <vector>

namespace ana {

namespace NDLArGeo {

// define the full TPC Dimensions from YAML file
static constexpr double TPC_DX = 46.788;
static constexpr double TPC_DY = 297.6;
static constexpr double TPC_DZ = 95.232;

static constexpr double HALF_DX = TPC_DX / 2.0;
static constexpr double HALF_DY = TPC_DY / 2.0;
static constexpr double HALF_DZ = TPC_DZ / 2.0;

// new definitions for detector geometry,
// using maximum boundaries of TPC modules
static constexpr double NDLArXLo = -276.2885 - HALF_DX;
static constexpr double NDLArXHi = 323.7115 + HALF_DX;
static constexpr double NDLArYLo = -66.8713 - HALF_DY;
static constexpr double NDLArYHi = -66.8713 + HALF_DY;
static constexpr double NDLArZLo = 465.7559 - HALF_DZ;
static constexpr double NDLArZHi = 865.7559 + HALF_DZ;

// Center positions of all 70 TPCs [X, Y, Z]
static const std::vector<std::vector<double>> TPC_CENTERS = {
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
    {-76.2885, -66.8713, 465.7559},  {-123.7115, -66.8713, 465.7559},
    {-76.2885, -66.8713, 565.7559},  {-123.7115, -66.8713, 565.7559},
    {-76.2885, -66.8713, 665.7559},  {-123.7115, -66.8713, 665.7559},
    {-76.2885, -66.8713, 765.7559},  {-123.7115, -66.8713, 765.7559},
    {-76.2885, -66.8713, 865.7559},  {-123.7115, -66.8713, 865.7559},
    {23.7115, -66.8713, 465.7559},   {-23.7115, -66.8713, 465.7559},
    {23.7115, -66.8713, 565.7559},   {-23.7115, -66.8713, 565.7559},
    {23.7115, -66.8713, 665.7559},   {-23.7115, -66.8713, 665.7559},
    {23.7115, -66.8713, 765.7559},   {-23.7115, -66.8713, 765.7559},
    {23.7115, -66.8713, 865.7559},   {-23.7115, -66.8713, 865.7559},
    {123.7115, -66.8713, 465.7559},  {76.2885, -66.8713, 465.7559},
    {123.7115, -66.8713, 565.7559},  {76.2885, -66.8713, 565.7559},
    {123.7115, -66.8713, 665.7559},  {76.2885, -66.8713, 665.7559},
    {123.7115, -66.8713, 765.7559},  {76.2885, -66.8713, 765.7559},
    {123.7115, -66.8713, 865.7559},  {76.2885, -66.8713, 865.7559},
    {223.7115, -66.8713, 465.7559},  {176.2885, -66.8713, 465.7559},
    {223.7115, -66.8713, 565.7559},  {176.2885, -66.8713, 565.7559},
    {223.7115, -66.8713, 665.7559},  {176.2885, -66.8713, 665.7559},
    {223.7115, -66.8713, 765.7559},  {176.2885, -66.8713, 765.7559},
    {223.7115, -66.8713, 865.7559},  {176.2885, -66.8713, 865.7559},
    {323.7115, -66.8713, 465.7559},  {276.2885, -66.8713, 465.7559},
    {323.7115, -66.8713, 565.7559},  {276.2885, -66.8713, 565.7559},
    {323.7115, -66.8713, 665.7559},  {276.2885, -66.8713, 665.7559},
    {323.7115, -66.8713, 765.7559},  {276.2885, -66.8713, 765.7559},
    {323.7115, -66.8713, 865.7559},  {276.2885, -66.8713, 865.7559}};

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

////////////////////////////////////////////////////////////////////////////////
// Position Checking functions:

// Checks if position (x, y, z) is inside the outer active volume box
inline bool
IsInOuterFiducialVolume(double x, double y, double z,
                        double margin_x_high = 0.0, double margin_x_low = 0.0,
                        double margin_y_high = 0.0, double margin_y_low = 0.0,
                        double margin_z_high = 0.0, double margin_z_low = 0.0) {

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

inline bool IsInCombinedNDLArFV(double x, double y, double z, double o_mx_high,
                                double o_mx_low, double o_my_high,
                                double o_my_low, double o_mz_high,
                                double o_mz_low, double box_margin) {

  if (!IsInOuterFiducialVolume(x, y, z, o_mx_high, o_mx_low, o_my_high,
                               o_my_low, o_mz_high, o_mz_low))
    return false;

  return IsInBoxesFiducialVolume(x, y, z, box_margin, box_margin, box_margin);
}

inline bool BoxContainedInBox(double x_min, double y_min, double z_min,
                              double x_max, double y_max, double z_max,
                              double Containment_crop) {
  return (IsInOuterFiducialVolume(x_min, y_min, z_min, Containment_crop,
                                  Containment_crop, Containment_crop,
                                  Containment_crop, Containment_crop,
                                  Containment_crop) &&
          IsInOuterFiducialVolume(x_max, y_max, z_max, Containment_crop,
                                  Containment_crop, Containment_crop,
                                  Containment_crop, Containment_crop,
                                  Containment_crop));
}

} // namespace NDLArGeo
} // namespace ana
#endif // NDLARGEOMETRY_H
