#include "duneanaobj/StandardRecord/SRInteraction.h"

namespace ana {

constexpr static auto kCannotFindTruth = std::numeric_limits<size_t>::max();

inline size_t GetTrueInteractionIndex(caf::SRInteraction const &nu_int) {
  if (!nu_int.truthOverlap.size()) {
    return kCannotFindTruth;
  }
  return std::distance(
      nu_int.truthOverlap.begin(),
      std::max_element(nu_int.truthOverlap.begin(), nu_int.truthOverlap.end()));
}
} // namespace ana
