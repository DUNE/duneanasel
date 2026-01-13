
if(duneanasel_USE_SRProxy)
  add_executable(separate_FD_samples separate_FD_samples.cxx)
  target_link_libraries(separate_FD_samples duneanasel::all)
  install(TARGETS
    separate_FD_samples
    DESTINATION bin)
endif()

add_executable(NDLAr_NumuCCContained_LongestTrack NDLAr_NumuCCContained_LongestTrack.cxx)
target_link_libraries(NDLAr_NumuCCContained_LongestTrack duneanasel::all)

install(TARGETS
  NDLAr_NumuCCContained_LongestTrack
  DESTINATION bin)
