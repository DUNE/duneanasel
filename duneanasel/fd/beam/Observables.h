#pragma once

#include "duneanasel/fd/beam/Selections.h"

#include <limits>

namespace proj::beam {

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
auto det_pos_cm(T const &fd_int) {
  return fd_int.vtx;
}

namespace fd1x8x6 {

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline float ENuReco(T const &fd_int, sel::beam::fd1x8x6::Sample smpl) {
  switch (smpl) {
  case sel::beam::fd1x8x6::kNuMuCCLike:
  case sel::beam::fd1x8x6::kNuECCLike:
    return fd_int.Enu.lep_calo;
  case sel::beam::fd1x8x6::kNCLike:
    return fd_int.Enu.calo;
  case sel::beam::fd1x8x6::kRejected:
  default:
    return std::numeric_limits<float>::quiet_NaN();
  }
}

} // namespace fd1x8x6

namespace FDVD = fd1x8x6;
namespace FD1 = FDVD;

} // namespace proj::beam
