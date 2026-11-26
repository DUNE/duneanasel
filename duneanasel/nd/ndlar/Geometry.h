#pragma once

#include "Math/Point3D.h"

#include <vector>

namespace geo {

namespace ndlar {

// LP: Add provenance info for these (which YAML file?)
const double TPC_DX = 46.788;
const double TPC_DY = 297.6;
const double TPC_DZ = 95.232;

const double HALF_DX = TPC_DX / 2.0;
const double HALF_DY = TPC_DY / 2.0;
const double HALF_DZ = TPC_DZ / 2.0;

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

  double DistanceToWall(double const &x, double const &y, double const &z) {
    return std::min({
        x - min.x(),
        max.x() - x,
        y - min.y(),
        max.y() - y,
        z - min.z(),
        max.z() - z,
    });
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

/// %brief Checks if a point is inside the cropped active volume box
inline bool PointInNDLArCroppedActiveVolume(
    double x, double y, double z, double x_high_crop = 0.0,
    double x_low_crop = 0.0, double y_high_crop = 0.0, double y_low_crop = 0.0,
    double z_high_crop = 0.0, double z_low_crop = 0.0) {

  static auto const NDLArAABB = GetNDLArAABB();

  return NDLArAABB.CroppedContains(x, y, z,
                                   AABBCrop{x_low_crop, x_high_crop, y_low_crop,
                                            y_high_crop, z_low_crop,
                                            z_high_crop});
}

/// %brief Checks if a point is inside any TPC volume with crops
inline bool PointInTPCCroppedVolume(double x, double y, double z,
                                    double x_crop = 0.0, double y_crop = 0.0,
                                    double z_crop = 0.0) {
  static auto const TPCAABBs = GetTPCAABBs();

  for (const auto &tpc : TPCAABBs) {
    if (tpc.CroppedContains(
            x, y, z,
            AABBCrop{x_crop, x_crop, y_crop, y_crop, z_crop, z_crop})) {
      return true;
    }
  }
  return false;
}

/// %brief Checks if a point is inside the union of the cropped active and
/// individual TPC volumes
inline bool PointInCroppedLArVolume(double x, double y, double z,
                                    double x_high_crop, double x_low_crop,
                                    double y_high_crop, double y_low_crop,
                                    double z_high_crop, double z_low_crop,
                                    double TPC_crop) {

  return PointInNDLArCroppedActiveVolume(x, y, z, x_high_crop, x_low_crop,
                                         y_high_crop, y_low_crop, z_high_crop,
                                         z_low_crop) &&
         PointInTPCCroppedVolume(x, y, z, TPC_crop, TPC_crop, TPC_crop);
}

} // namespace ndlar
} // namespace geo
