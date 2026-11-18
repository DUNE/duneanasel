#ifndef PHYSICAL_CUTS
#define PHYSICAL_CUTS

// reco
bool HasSandParticles(const caf::SRInteractionProxy* recoint)
{
  return recoint->part.sandreco.size() > 0;
}

int NPrimaryMuons(const caf::SRInteractionProxy* recoint)
{
  int n=0;
  for (const auto& recopart : recoint->part.sandreco) {
    if (recopart.primary == 1 && std::abs(recopart.pdg) == 13) {++n;}
  }
  return n;
}

#endif

