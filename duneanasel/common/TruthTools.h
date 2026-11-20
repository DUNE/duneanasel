#include "duneanaobj/StandardRecord/SRInteraction.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

#include "duneanasel/common/type_traits.h"

namespace ana {

constexpr static auto kCannotFindTruth = std::numeric_limits<size_t>::max();

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline size_t GetTrueInteractionIndex(T const &nu_int) {
  if (!nu_int.truthOverlap.size()) {
    return kCannotFindTruth;
  }

  size_t idx = 0;
  double max_el = std::numeric_limits<double>::lowest();
  for(int i = 0; i < nu_int.truthOverlap.size(); ++i){
    if(nu_int.truthOverlap[i] > max_el){
      max_el = nu_int.truthOverlap[i];
      idx = i;
    }
  }
  if (idx >= nu_int.truth.size()) {
    return kCannotFindTruth;
  }
  return nu_int.truth[idx];
}
} // namespace ana
