#include "SRInteraction.h"

namespace sel {
namespace atm {

enum Sample { kRejected = 0, kNuMuCCLike, kNuECCLike, kNCLike };

namespace fd1x2x6 {

// private correspondance from Pierre Granger
inline bool InFV(caf::SRInteraction const &fd_int) {
  if (std::fabs(fd_int.vtx.x) > 343) {
    return false;
  }
  if (std::fabs(fd_int.vtx.y) > 588) {
    return false;
  }
  float const &z = fd_int.vtx.z;
  if ((z < 20) || (z > 1374)) {
    return false;
  }

  return true;
}

// https://github.com/DUNE/atmospherics-tools/blob/401037e4ced8ff1ab48fefacf0929bbc260b6489/weight_calculation/app/split_channels.cxx#L142-L153
inline bool NumuCCLike(caf::SRInteraction const &fd_int) {
  return fd_int.nuhyp.cvn.numu > 0.56;
}

inline bool NueCCLike(caf::SRInteraction const &fd_int) {

  return !NumuCCLike(fd_int) && (fd_int.nuhyp.cvn.nue > 0.56);
}

inline bool NCLike(caf::SRInteraction const &fd_int) {
  return !NumuCCLike(fd_int) && !NueCCLike(fd_int);
}

inline Sample ApplySelection(caf::SRInteraction const &fd_int) {
  if (sel::atm::fd1x2x6::NuMuCCLike(fd_int)) {
    return sel::atm::kNuMuCCLike;
  } else if (sel::atm::fd1x2x6::NuMuCCLike(fd_int)) {
    return sel::atm::kNuECCLike;
  } else if (sel::atm::fd1x2x6::NCLike(fd_int)) {
    return sel::atm::kNCLike;
  }
  return sel::atm::kRejected;
}

} // namespace fd1x2x6
} // namespace atm
} // namespace sel
