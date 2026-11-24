#pragma once

#include <algorithm>
#include <cmath>

#include "duneanasel/common/RecoTools.h"
#include "duneanasel/nd/ndlar/NDLArGeometry.h"

// The selection in this file was originally developed by Eva Sabater; what
// follows is a transcription from the original CAFAna code.
//
// This numu CC inclusive selection header extends Selections.h
// Uses NDLArGeometry.h adapted from Freddi Vosper's version

namespace sel{
    namespace beam{

        namespace ndlar{

            enum Sample {
                kRejected = 0,
                kNuMuCCLikeContained,
                kNuECCLikeContained,
                kNuMuCCInclusive
            };

            // previous ND-LAr geometry definitions
            // AV outer boundaries [cm]
            constexpr double NDLArXLo = -346.9;
            constexpr double NDLArXHi = 346.9;
            constexpr double NDLArYLo = -215.5;
            constexpr double NDLArYHi = 81.7;
            constexpr double NDLArZLo = 418.2;
            constexpr double NDLArZHi = 913.3;

            // functions form Selections.h
            template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
            inline bool IntHasParticles(T const &nd_int) {
                if (nd_int.part.dlp.empty()) return false;
                return true;
            }

            template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
            inline bool LongestPrimTrackLengthCut(T const &nd_int) {
                auto plong = ana::GetLongestParticle(nd_int);
                if (!plong) return false;
                return ana::ParticleLength(*plong) > 0.1;
            }

            template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
            inline bool InFV(T const &nd_int) {
                float x = nd_int.vtx.x;
                if ((x <= (NDLArXLo + 25)) || (x >= (NDLArXHi - 25))) return false;
                float y = nd_int.vtx.y;
                if ((y <= (NDLArYLo + 25)) || (y >= (NDLArYHi - 25))) return false;
                float z = nd_int.vtx.z;
                if ((z <= (NDLArZLo + 25)) || (z >= (NDLArZHi - 25))) return false;
                return true;
            }

            template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
            inline bool AllPrimaryParticlesContained(T const &nd_int) {
                for (auto const &p : nd_int.part.dlp)
                    if ((p.primary == 1) && (p.contained != 1)) return false;
                return true;
            }

            template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
            inline bool ParticlesNDLArContainedMuonsEscapeDownstream(T const &nd_int) {
                for (auto const &p : nd_int.part.dlp) {
                    bool contained    = (p.contained == 1);
                    bool isPrimaryMuon = (p.pdg == 13 && p.primary == 1);
                    bool escapesToTMS  = false;
                    if (isPrimaryMuon) {
                    escapesToTMS =
                        (p.end.x < (NDLArXHi - 25)) && (p.end.x > (NDLArXLo + 25)) &&
                        (p.end.y < (NDLArYHi - 25)) && (p.end.y > (NDLArYLo + 25)) &&
                        (p.end.z > (NDLArZLo + 25));
                    }
                    if (!contained && !escapesToTMS) return false;
                }
                return true;
            }

            template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
            inline bool HasParticleWithReconstructedPID(T const &nd_int, int pid) {
                for (auto const &p : nd_int.part.dlp)
                    if ((p.primary == 1) && (p.pdg == pid)) return true;
                return false;
            }

            namespace numode{
                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool NuMuCCLikeContained(T const &nd_int) {
                    if (!InFV(nd_int) || !AllPrimaryParticlesContained(nd_int)) return false;
                    auto plong = ana::GetLongestParticle(nd_int);
                    if (!plong || (plong->pdg != 13)) return false;
                    return true;
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool NuMuCCLike(T const &nd_int) {
                    for (auto const &p : nd_int.part.dlp)
                        if (p.primary == 1 && p.pdg == 13) return true;
                    return false;
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool NuECCLikeContained(T const &nd_int) {
                    if (!InFV(nd_int) || !AllPrimaryParticlesContained(nd_int) ||
                        !HasParticleWithReconstructedPID(nd_int, 11)) return false;
                    return true;
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline sel::beam::ndlar::Sample ApplySelectionV0(T const &fd_int) {
                    if (!sel::beam::ndlar::InFV(fd_int)) return sel::beam::ndlar::kRejected;
                    if (sel::beam::ndlar::numode::NuMuCCLikeContained(fd_int))
                        return sel::beam::ndlar::kNuMuCCLikeContained;
                    else if (sel::beam::ndlar::numode::NuECCLikeContained(fd_int))
                        return sel::beam::ndlar::kNuECCLikeContained;
                    return sel::beam::ndlar::kRejected;
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline sel::beam::ndlar::Sample ApplySelectionV0p1(T const &nd_int) {
                    if (!sel::beam::ndlar::IntHasParticles(nd_int) ||
                        !sel::beam::ndlar::LongestPrimTrackLengthCut(nd_int) ||
                        !sel::beam::ndlar::InFV(nd_int) ||
                        !sel::beam::ndlar::ParticlesNDLArContainedMuonsEscapeDownstream(nd_int))
                        return sel::beam::ndlar::kRejected;
                    if (sel::beam::ndlar::numode::NuMuCCLike(nd_int))
                        return sel::beam::ndlar::kNuMuCCLikeContained;
                    else if (sel::beam::ndlar::numode::NuECCLikeContained(nd_int))
                        return sel::beam::ndlar::kNuECCLikeContained;
                    return sel::beam::ndlar::kRejected;
                }

            } // numode

            // numu CC inclusive selection
            namespace numuCC{

                // gemoetrical boundaries in [cm]
                // values come from FV study done by Freddie
                constexpr double kFV_OX_Hi = 10.;
                constexpr double kFV_OX_Lo = 10.;
                constexpr double kFV_OY_Hi = 10.;
                constexpr double kFV_OY_Lo = 10.;
                constexpr double kFV_OZ_Hi = 20.;
                constexpr double kFV_OZ_Lo = 40.;
                constexpr double kFV_Box =  2.;
                constexpr double kMuonMinLen = 25.;
                constexpr double kGeomContMargin =  5.;
                constexpr double kWallMargin =  5.;

                template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
                inline double ParticleLen(T const &p) {
                    if (std::isnan(p.start.x) || std::isinf(p.start.x) ||
                        std::isnan(p.end.x)   || std::isinf(p.end.x))   return -1.;
                    if (std::isnan(p.start.y) || std::isinf(p.start.y) ||
                        std::isnan(p.end.y)   || std::isinf(p.end.y))   return -1.;
                    if (std::isnan(p.start.z) || std::isinf(p.start.z) ||
                        std::isnan(p.end.z)   || std::isinf(p.end.z))   return -1.;
                    return std::hypot(p.end.x - p.start.x,
                                    p.end.y - p.start.y,
                                    p.end.z - p.start.z);
                }

                inline double DistToOuterWall(double x, double y, double z) {
                    double dx = std::min(x - geo::kNDLArXLo, geo::kNDLArXHi - x);
                    double dy = std::min(y - geo::kNDLArYLo, geo::kNDLArYHi - y);
                    double dz = std::min(z - geo::kNDLArZLo, geo::kNDLArZHi - z);
                    return std::min({dx, dy, dz});
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool VtxInCombinedFV(T const &nd_int) {
                    return geo::IsInCombinedFV(nd_int.vtx.x, nd_int.vtx.y, nd_int.vtx.z,
                                                kFV_OX_Hi, kFV_OX_Lo,
                                                kFV_OY_Hi, kFV_OY_Lo,
                                                kFV_OZ_Hi, kFV_OZ_Lo,
                                                kFV_Box);
                }

                template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
                inline bool IsPrimaryMuon(T const &p) {
                    return (p.primary == 1) && (p.pdg == 13 || p.pdg == -13);
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool ExactlyOnePrimaryMuon(T const &nd_int) {
                    int n = 0;
                    for (auto const &p : nd_int.part.dlp)
                        if (IsPrimaryMuon(p)) ++n;
                    return n == 1;
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool MuonLengthCut(T const &nd_int) {
                    for (auto const &p : nd_int.part.dlp)
                        if (IsPrimaryMuon(p) && ParticleLen(p) >= kMuonMinLen) return true;
                    return false;
                }

                template <typename T, typename C = Proxyable_t<caf::SRRecoParticle, T>>
                inline bool MuonEscapesDownstream(T const &p) {
                    if (!(p.start.x > geo::kNDLArXLo + kWallMargin &&
                        p.start.x < geo::kNDLArXHi - kWallMargin &&
                        p.start.y > geo::kNDLArYLo + kWallMargin &&
                        p.start.y < geo::kNDLArYHi - kWallMargin &&
                        p.start.z > geo::kNDLArZLo + kWallMargin &&
                        p.start.z < geo::kNDLArZHi - kWallMargin)) return false;
                    return (p.end.x > geo::kNDLArXLo && p.end.x < geo::kNDLArXHi &&
                            p.end.y > geo::kNDLArYLo && p.end.y < geo::kNDLArYHi &&
                            p.end.z > geo::kNDLArZLo);
                }

                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool SPINEContained(T const &nd_int) {
                    for (auto const &p : nd_int.part.dlp) {
                        if (p.primary != 1) continue;
                        if (IsPrimaryMuon(p)) {
                            if (p.contained != 1 && !MuonEscapesDownstream(p)) return false;
                        } else {
                            if (p.contained != 1) return false;
                        }
                    }
                    return true;
                }

                // preliminary geomatrical containment
                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline bool GeomContained(T const &nd_int) {
                    for (auto const &p : nd_int.part.dlp) {
                        if (p.primary != 1) continue;
                        if (IsPrimaryMuon(p)) {
                            if (!MuonEscapesDownstream(p)) return false;
                        } else {
                            if (std::isnan(p.end.x) || std::isinf(p.end.x)) return false;
                            if (DistToOuterWall(p.end.x, p.end.y, p.end.z) < kGeomContMargin)
                                return false;
                        }
                    }
                    return true;
                }

                // selection using SPINE containment
                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline sel::beam::ndlar::Sample ApplySelectionV1_SPINE(T const &nd_int) {
                    if (!VtxInCombinedFV(nd_int))       return sel::beam::ndlar::kRejected;
                    if (!ExactlyOnePrimaryMuon(nd_int)) return sel::beam::ndlar::kRejected;
                    if (!MuonLengthCut(nd_int))         return sel::beam::ndlar::kRejected;
                    if (!SPINEContained(nd_int))        return sel::beam::ndlar::kRejected;
                    return sel::beam::ndlar::kNuMuCCInclusive;
                }

                // selection using distance to wall (geometry) containment
                template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
                inline sel::beam::ndlar::Sample ApplySelectionV1_Geom(T const &nd_int) {
                    if (!VtxInCombinedFV(nd_int))       return sel::beam::ndlar::kRejected;
                    if (!ExactlyOnePrimaryMuon(nd_int)) return sel::beam::ndlar::kRejected;
                    if (!MuonLengthCut(nd_int))         return sel::beam::ndlar::kRejected;
                    if (!GeomContained(nd_int))         return sel::beam::ndlar::kRejected;
                    return sel::beam::ndlar::kNuMuCCInclusive;
                }
            } //numuCC
        } // ndlar
    } // beam
} // sel
