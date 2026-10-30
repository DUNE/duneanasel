add_executable(ndlarccnumu0pievl ndlarccnumu0pievl.cc)
target_link_libraries(ndlarccnumu0pievl duneanasel::all duneanaobj::StandardRecordProxy)

install(TARGETS ndlarccnumu0pievl DESTINATION bin)
