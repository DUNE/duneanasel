add_executable(Event_Loop Event_Loop.cc)
target_link_libraries(Event_Loop duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

# add_executable(Nu_mu_CC0Pi_RNTuple_read Nu_mu_CC0Pi_RNTuple_read.cc)
# target_link_libraries(Nu_mu_CC0Pi_RNTuple_read duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

install(TARGETS
Event_Loop
# Nu_mu_CC0Pi_RNTuple_read
DESTINATION bin)
