#include <Math/Vector3D.h>
#include <ROOT/RVec.hxx>

using V3 = ROOT::Math::XYZVector;
template <typename T>
using RVec = ROOT::VecOps::RVec<T>;

inline RVec<double> TrackLengths(const RVec<RVec<V3>> &starts,
                                const RVec<RVec<V3>> &ends) {
  RVec<double> lengths;
  for (std::size_t i = 0; i < starts.size(); ++i) {
    for (std::size_t j = 0; j < starts[i].size(); ++j) {
      lengths.push_back((ends[i][j] - starts[i][j]).R());
    }
  }
  return lengths;
}