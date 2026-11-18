#ifndef GEOMETRY_H
#define GEOMETRY_H

const double kSandOffsetY = 238.473;   // cm
const double kSandOffsetZ = 2391.0;    // cm

inline caf::SRVector3D NDHallToSANDCoordinates(const caf::SRVector3DProxy& ndhall_coords){
  // Convert coordinates from ND Hall to SAND coordinates
  // z points along beam, y points up, x points to the left if you look at the beam direction
  caf::SRVector3D sand_coords;
  sand_coords.x = ndhall_coords.x;                      // x remains the same
  sand_coords.y = ndhall_coords.y + kSandOffsetY;       // cm   
  sand_coords.z = ndhall_coords.z - kSandOffsetZ;       // cm   

  return sand_coords;
}

inline caf::SRVector3D SANDToNDHallCoordinates(const caf::SRVector3DProxy& sand_coords){
  // Convert coordinates from SAND to ND Hall coordinates
  // z points along beam, y points up, x points to the left if you look at the beam direction
  caf::SRVector3D ndhall_coords;
  ndhall_coords.x = sand_coords.x;                      // x remains the same
  ndhall_coords.y = sand_coords.y - kSandOffsetY;       // cm   
  ndhall_coords.z = sand_coords.z + kSandOffsetZ;       // cm
  return ndhall_coords;
}

bool inFrontMagnet(const caf::SRInteractionProxy* recoint){
  auto nu_vtx = NDHallToSANDCoordinates(recoint->vtx);

  // auto radius = sqrt(std::pow(nu_vtx.y, 2) + std::pow(nu_vtx.z, 2));

  if (std::abs(nu_vtx.x) < 169 && std::abs(nu_vtx.y) < 200) {
    return true;
  }

  return false;
}

// fiducial cut for the front ECAL:
// considering the center of SAND, we get a region with 2.0 < r < 2.3 m
// and z < 0.910 (which actually is z < 23 m in global NDHall coordinates)
bool inFrontECAL(const caf::SRInteractionProxy* recoint){
  auto nu_vtx = NDHallToSANDCoordinates(recoint->vtx);
  auto radius = sqrt(std::pow(nu_vtx.y, 2) + std::pow(nu_vtx.z, 2));

  auto z_cut_glob = 2300.0; // cm
  auto z_cut_loc = z_cut_glob - 2391.0;

  if (radius < 230 && radius > 200 && nu_vtx.z < z_cut_loc && std::abs(nu_vtx.x) < 169) {
    return true;
  }

  return false;
}

bool inSANDInnerVol(const caf::SRInteractionProxy* recoint){
  auto nu_vtx = NDHallToSANDCoordinates(recoint->vtx);
  auto radius = sqrt(std::pow(nu_vtx.y, 2) + std::pow(nu_vtx.z, 2));

  if (radius < 190 && std::abs(nu_vtx.x) < 159) {
    return true;
  }

  return false;
}

bool inSANDFV(const caf::SRInteractionProxy* recoint){
  return inSANDInnerVol(recoint) || inFrontECAL(recoint);
}

#endif