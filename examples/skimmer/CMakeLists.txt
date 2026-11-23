add_executable(skimmer skimmer.cxx)
target_link_libraries(skimmer
  duneanasel::all)

add_executable(skimmer_rdf skimmer_rdf.cxx)
target_link_libraries(skimmer_rdf
  duneanasel::all
  ROOT::Core
  ROOT::GenVector
  ROOT::ROOTVecOps
  ROOT::ROOTDataFrame
  ROOT::RIO
  ROOT::Tree
  ROOT::ROOTNTuple
  ROOT::Hist
  ROOT::Gpad
  ROOT::Graf
)

add_executable(skimmer_rntuple skimmer_rntuple.cxx)
target_link_libraries(skimmer_rntuple
  duneanasel::all
  ROOT::Core
  ROOT::GenVector
  ROOT::ROOTVecOps
  ROOT::ROOTDataFrame
  ROOT::RIO
  ROOT::Tree
  ROOT::ROOTNTuple
  ROOT::Hist
  ROOT::Gpad
  ROOT::Graf
)