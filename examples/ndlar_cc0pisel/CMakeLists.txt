add_executable(Event_Loop Event_Loop.cc)
target_link_libraries(Event_Loop duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Nu_mu_CC0Pi_RNTuple_read Nu_mu_CC0Pi_RNTuple_read.cc)
target_link_libraries(Nu_mu_CC0Pi_RNTuple_read duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Nu_mu_CC0Pi_RNTuple_opt Nu_mu_CC0Pi_RNTuple_opt.cc)
target_link_libraries(Nu_mu_CC0Pi_RNTuple_opt duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Opt_param_plots Opt_param_plots.cc)
target_link_libraries(Opt_param_plots duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Stacked_plot_check Stacked_plot_check.cc)
target_link_libraries(Stacked_plot_check duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(PlotFixedCuts PlotFixedCuts.cc)
target_link_libraries(PlotFixedCuts duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(PlotFixedCuts_andFV PlotFixedCuts_andFV.cc)
target_link_libraries(PlotFixedCuts_andFV duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Plot_Reco Plot_Reco.cc)
target_link_libraries(Plot_Reco duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Plot_True Plot_True.cc)
target_link_libraries(Plot_True duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Plot_edges Plot_edges.cc)
target_link_libraries(Plot_edges duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Reco_with_muon_mom Reco_with_muon_mom.cc)
target_link_libraries(Reco_with_muon_mom duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)

add_executable(Nu_mu_CC0Pi_Apply_Sel Nu_mu_CC0Pi_Apply_Sel.cc)
target_link_libraries(Nu_mu_CC0Pi_Apply_Sel duneanasel::all duneanaobj::StandardRecordProxy ROOT::ROOTNTuple)


install(TARGETS
Event_Loop
Nu_mu_CC0Pi_RNTuple_read
Nu_mu_CC0Pi_RNTuple_opt
Opt_param_plots
Stacked_plot_check
PlotFixedCuts
PlotFixedCuts_andFV
Plot_Reco
Plot_True
Plot_edges
Reco_with_muon_mom
Nu_mu_CC0Pi_Apply_Sel 
DESTINATION bin)
