#pragma once

#include "duneanaobj/StandardRecord/SRInteraction.h"

#include "duneanasel/common/type_traits.hxx"

namespace sel {
namespace beam {

enum Sample { kRejected = 0, kNuMuCCLike, kNuECCLike, kNCLike };

namespace fd1x8x6 {

// from CVN training, here:
// https://github.com/DUNE/dunereco/blob/ff92d1bc966c6e3cf7db1a434c8a723e74b0712c/dunereco/CVN/art/CVNZlibMaker_module.cc#L295
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool InFV(T const &fd_int) {

  if (std::fabs(fd_int.vtx.x) > 300) {
    return false;
  }
  if (std::fabs(fd_int.vtx.y) > 680) {
    return false;
  }
  float z = fd_int.vtx.z;
  if ((z < 50) || (z > 850)) {
    return false;
  }

  return true;
}

namespace numode {

// Motivating analysis in jupyter notebook in this repository @ /examples/cut_optimization/fd/VDFDToySelections.ipynb
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool NuMuCCLike(T const &fd_int) {

  return (fd_int.nuhyp.cvn.numu > 0.21) && (fd_int.nuhyp.cvn.nue < 0.77) &&
         (fd_int.nuhyp.cvn.nc < 0.54);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool NuECCLike(T const &fd_int) {

  return (fd_int.nuhyp.cvn.nue > 0.58) && (fd_int.nuhyp.cvn.numu < 0.051) &&
         (fd_int.nuhyp.cvn.nc < 0.071);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool NCLike(T const &fd_int) {

  return (fd_int.nuhyp.cvn.nc > 0.32) && (fd_int.nuhyp.cvn.numu < 0.34) &&
         (fd_int.nuhyp.cvn.nue < 0.68);
}

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline sel::beam::Sample ApplySelection(T const &fd_int) {

  if (sel::beam::fd1x8x6::numode::NuMuCCLike(fd_int)) {
    return sel::beam::kNuMuCCLike;
  } else if (sel::beam::fd1x8x6::numode::NuMuCCLike(fd_int)) {
    return sel::beam::kNuECCLike;
  } else if (sel::beam::fd1x8x6::numode::NCLike(fd_int)) {
    return sel::beam::kNCLike;
  }
  return sel::beam::kRejected;
}

} // namespace numode
} // namespace fd1x8x6

namespace fd1x2x6 {
// from CVN training, here:
//
// https://github.com/DUNE/dunereco/blob/ff92d1bc966c6e3cf7db1a434c8a723e74b0712c/dunereco/CVN/art/CVNZlibMaker_module.cc#L295
template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool InFV(T const &fd_int) {

  if (std::fabs(fd_int.vtx.X()) > 310) {
    return false;
  }
  if (std::fabs(fd_int.vtx.Y()) > 550) {
    return false;
  }
  auto const &z = fd_int.vtx.Z();
  if ((z < 50) || (z > 1244)) {
    return false;
  }

  return true;
}

} // namespace fd1x2x6

} // namespace beam
} // namespace sel
