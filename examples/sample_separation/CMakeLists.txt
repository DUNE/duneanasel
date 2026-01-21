if(duneanasel_USE_SRProxy)
  add_executable(separate_FD_samples separate_FD_samples.cxx)
  target_link_libraries(separate_FD_samples duneanasel::all)
  install(TARGETS
    separate_FD_samples
    DESTINATION bin)
endif()
