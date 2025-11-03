#include "duneanasel/fd/beam/Selections.h"

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"

#include <memory>

int main(int argc, char const *argv[]) {

  TChain pch("cafmaker/cafTree");
  pch.Add(argv[1]);
  caf::StandardRecordProxy srp(&pch, "rec");

  Long64_t ents = pch.GetEntries();
  std::cout << "Input tree has " << ents << " entries." << std::endl;

  for (Long64_t i = 0; i < ents; ++i) {
    pch.GetEntry(i);

    //analysis goes here

  }
}
