add_executable(NDLAr_NumuCCContained_LongestTrack NDLAr_NumuCCContained_LongestTrack.cxx)
target_link_libraries(NDLAr_NumuCCContained_LongestTrack
  duneanasel::all duneanaobj::StandardRecordProxy)

add_executable(NDLAr_NumuCCContained_VtxPosition NDLAr_NumuCCContained_VtxPosition.cxx)
target_link_libraries(NDLAr_NumuCCContained_VtxPosition
  duneanasel::all duneanaobj::StandardRecordProxy)

add_executable(track_lengths_sr track_lengths_sr.cxx)
target_link_libraries(track_lengths_sr
  duneanasel::all duneanaobj::StandardRecordProxy)

add_executable(track_lengths_srproxy track_lengths_srproxy.cxx)
target_link_libraries(track_lengths_srproxy
  duneanasel::all duneanaobj::StandardRecordProxy)

add_executable(track_lengths_rdf track_lengths_rdf.cxx)
target_link_libraries(track_lengths_rdf
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

add_executable(track_lengths_rdf_dumb track_lengths_rdf_dumb.cxx)
target_link_libraries(track_lengths_rdf_dumb
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

add_executable(track_lengths_rdf_super_dumb track_lengths_rdf_super_dumb.cxx)
target_link_libraries(track_lengths_rdf_super_dumb
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

add_executable(track_lengths_rdf_super_super_dumb track_lengths_rdf_super_super_dumb.cxx)
target_link_libraries(track_lengths_rdf_super_super_dumb
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

install(TARGETS
  NDLAr_NumuCCContained_LongestTrack
  NDLAr_NumuCCContained_VtxPosition
  track_lengths_sr
  track_lengths_srproxy
  track_lengths_rdf
  track_lengths_rdf_dumb
  track_lengths_rdf_super_dumb
  track_lengths_rdf_super_super_dumb
  DESTINATION bin)
