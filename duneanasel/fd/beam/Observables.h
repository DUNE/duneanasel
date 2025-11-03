#pragma once

#include "duneanasel/fd/beam/Selections.h"

#include <limits>

namespace proj::beam::fd1x8x6 {

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline float ENuReco(T const &fd_int, sel::beam::Sample smpl) {
  switch (smpl) {
  case sel::kRejected:
    return std::numeric_limits<float>::quiet_NaN();
  case sel::kNuMuCCLike:
  case sel::kNuECCLike:
    return fd_int.Enu.lep_calo;
  case sel::kNCLike:
    return fd_int.Enu.calo;
  }
}

} // namespace proj::beam::fd1x8x6



