// ND-LAr geometry helper functions
// Adapted from Frddie Vosper's geometry header

#pragma once

#include <array>
#include <cmath>
#include <vector>

namespace sel{
    namespace beam{
        namespace ndlar{
            namespace geo{

                // TPC dimensions [cm] from yaml
                constexpr double TPC_DX = 46.788;
                constexpr double TPC_DY = 297.6;
                constexpr double TPC_DZ = 95.232;

                constexpr double HALF_DX = TPC_DX / 2.;
                constexpr double HALF_DY = TPC_DY / 2.;
                constexpr double HALF_DZ = TPC_DZ / 2.;

                // Outer boundaries [cm] (active volume)
                constexpr double NDLArXLo = -276.2885 - HALF_DX;
                constexpr double NDLArXHi = 323.7115 + HALF_DX;
                constexpr double NDLArYLo = -66.8713 - HALF_DY;
                constexpr double NDLArYHi = -66.8713 + HALF_DY;
                constexpr double NDLArZLo = 465.7559 - HALF_DZ;
                constexpr double NDLArZHi = 865.7559 + HALF_DZ;

                // Individual TPC box
                struct Box{
                    double x_min, x_max;
                    double y_min, y_max;
                    double z_min, z_max;

                    bool Contains(double x, double y, double z) const {
                        return (x >= x_min && x <= x_max &&
                                y >= y_min && y <= y_max &&
                                z >= z_min && z <= z_max);
                    }
                };

                // centre positions of all 70 TPCs [cm] {x, y, z}
                inline const std::vector<std::array<double,3>> TPCCentres = {
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
                    { 323.7115, -66.8713, 865.7559}, { 276.2885, -66.8713, 865.7559},
                };

                // TPC boxes with per-face margin
                inline std::vector<Box> GetTPCBoxes(double margin_x = 0., double margin_y = 0., double margin_z = 0.) {
                    std::vector<Box> boxes;
                    boxes.reserve(TPCCentres.size());
                    for (auto const &c : TPCCentres) {
                        boxes.push_back({c[0] - HALF_DX + margin_x,
                                         c[0] + HALF_DX - margin_x,
                                         c[1] - HALF_DY + margin_y,
                                         c[1] + HALF_DY - margin_y,
                                         c[2] - HALF_DZ + margin_z,
                                         c[2] + HALF_DZ - margin_z});
                    }
                    return boxes;
                }

                // can add a margin to outer active volume
                inline bool IsInOuterAV(double x, double y, double z,
                                        double mx_hi = 0., double mx_lo = 0.,
                                        double my_hi = 0., double my_lo = 0.,
                                        double mz_hi = 0., double mz_lo = 0.) {
                    return (x > NDLArXLo + mx_lo && x < NDLArXHi - mx_hi &&
                            y > NDLArYLo + my_lo && y < NDLArYHi - my_hi &&
                            z > NDLArZLo + mz_lo && z < NDLArZHi - mz_hi);
                }

                inline bool IsInAnyTPC(double x, double y, double z,
                                       double margin_x = 0., double margin_y = 0.,
                                       double margin_z =0.) {
                    for (auto const &box : GetTPCBoxes(margin_x, margin_y, margin_z)) {
                        if (box.Contains(x, y, z)) return true;
                    }
                    return false;
                }

                inline bool IsInCombinedFV(double x, double y, double z,
                                           double o_mx_hi, double o_mx_lo,
                                           double o_my_hi, double o_my_lo,
                                           double o_mz_hi, double o_mz_lo,
                                           double box_margin) {
                    if (!IsInOuterAV(x, y, z,
                                     o_mx_hi, o_mx_lo,
                                     o_my_hi, o_my_lo,
                                     o_mz_hi, o_mz_lo)) return false;
                    return IsInAnyTPC(x, y, z, box_margin, box_margin, box_margin);
                }

            } // geo
        } // ndlar
    } // beam
} // sel
