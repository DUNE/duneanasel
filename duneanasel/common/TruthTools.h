#include "duneanasel/common/type_traits.h"

#include <algorithm>
#include <optional>

namespace ana {

constexpr static auto kCannotFindTruth = std::numeric_limits<size_t>::max();

// template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
// inline size_t GetTrueInteractionIndex(T const &nu_int) {
//   if (!nu_int.truthOverlap.size()) {
//     return kCannotFindTruth;
//   }
//   size_t idx = std::distance(
//       nu_int.truthOverlap.begin(),
//       std::max_element(nu_int.truthOverlap.begin(), nu_int.truthOverlap.end()));
//   if (idx >= nu_int.truth.size()) {
//     return kCannotFindTruth;
//   }
//   return nu_int.truth[idx];
// }

// // uses Ancestor so only works for Proxy
// inline caf::SRTrueInteractionProxy *
// GetTrueInteraction(caf::SRInteractionProxy const &nu_int) {
//   if (nu_int.truth.empty()) {
//     return nullptr;
//   }

//   auto tidx = GetTrueInteractionIndex(nu_int);

//   if (tidx == kCannotFindTruth) {
//     return nullptr;
//   }

//   return &nu_int.Ancestor<caf::SRProxy>()->mc.nu[tidx];
// }

} // namespace ana
