#pragma once

#include "duneanasel/nd/ndlar/numucc/Constants.h"

namespace sel::beam::ndlar::numucc {

template <typename T, typename C = Proxyable_t<caf::SRInteraction, T>>
inline bool VtxInCombinedFV(T const &nd_int) {
  return geo::ndlar::PointInCroppedLArVolume(
      nd_int.vtx.x, nd_int.vtx.y, nd_int.vtx.z, kNDLAr_Crop_X_Hi_cm,
      kNDLAr_Crop_X_Lo_cm, kNDLAr_Crop_Y_Hi_cm, kNDLAr_Crop_Y_Lo_cm,
      kNDLAr_Crop_Z_Hi_cm, kNDLAr_Crop_Z_Lo_cm, kTPC_Crop_cm);
}

} // namespace sel::beam::ndlar::numucc
