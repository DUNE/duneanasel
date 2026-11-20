add_executable(MaCh3Winter26Skimmer MaCh3Winter26Skimmer.cc)
target_link_libraries(MaCh3Winter26Skimmer duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

install(TARGETS MaCh3Winter26Skimmer DESTINATION bin)
