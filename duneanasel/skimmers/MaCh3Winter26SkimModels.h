#pragma once

#include "Math/Point3D.h"
#include "ROOT/RNTupleModel.hxx"

#include <memory>
#include <string>
#include <vector>

namespace mach3winter26 {

struct Meta {
  std::unique_ptr<ROOT::RNTupleModel> model;

  std::shared_ptr<std::string> duneanaobj_input_version;
  std::shared_ptr<std::string> duneanaobj_process_version;
  std::shared_ptr<std::string> duneanasel_version;
  std::shared_ptr<double> pot_processed;
  std::shared_ptr<std::string> reco_path;
  std::shared_ptr<std::vector<std::string>> truth_categories;
  std::shared_ptr<std::vector<std::string>> truth_cctopologies;
  std::shared_ptr<std::vector<std::string>> param_names;
  std::shared_ptr<std::vector<std::vector<float>>> param_knots;

  Meta() {
    model = ROOT::RNTupleModel::Create();

    duneanaobj_input_version =
        model->MakeField<std::string>("duneanaobj_input_version");
    duneanaobj_process_version =
        model->MakeField<std::string>("duneanaobj_process_version");
    duneanasel_version = model->MakeField<std::string>("duneanasel_version");
    pot_processed = model->MakeField<double>("pot_processed");
    reco_path = model->MakeField<std::string>("reco_path");
    truth_categories =
        model->MakeField<std::vector<std::string>>("truth_categories");
    truth_cctopologies =
        model->MakeField<std::vector<std::string>>("truth_cctopologies");
    param_names = model->MakeField<std::vector<std::string>>("param_names");
    param_knots =
        model->MakeField<std::vector<std::vector<float>>>("param_knots");
  }
};

struct Event {
  std::unique_ptr<ROOT::RNTupleModel> model;

  std::shared_ptr<int> sr_entry;
  std::shared_ptr<int> ixn_entry;

  enum Sample { kRejected = 0, kNumuCCLike, kNueCCLike, kNCLike };
  struct {
    std::shared_ptr<int> sample;
    std::shared_ptr<ROOT::Math::XYZPoint> vtx_pos_cm;
    std::shared_ptr<double> enu_GeV;
    std::shared_ptr<double> mu_e_GeV;
    std::shared_ptr<double> shower_e_GeV;
    std::shared_ptr<double> sum_hadtrack_e_GeV;
  } reco;

  struct {
    std::shared_ptr<int> category;
    std::shared_ptr<int> cctopology;
    std::shared_ptr<bool> is_cc;
    std::shared_ptr<int> mode;
    std::shared_ptr<ROOT::Math::XYZPoint> vtx_pos_cm;
    std::shared_ptr<int> pdgnu;
    std::shared_ptr<int> pdgnu_original;
    std::shared_ptr<double> enu_GeV;
    std::shared_ptr<std::vector<std::vector<float>>> param_responses;
  } truth;

  Event() {

    model = ROOT::RNTupleModel::Create();

    sr_entry = model->MakeField<int>("sr_entry");
    ixn_entry = model->MakeField<int>("ixn_entry");

    reco.sample = model->MakeField<int>("sample");
    reco.vtx_pos_cm = model->MakeField<ROOT::Math::XYZPoint>("vtx_pos_cm");
    reco.enu_GeV = model->MakeField<double>("enu_rec_GeV");
    reco.mu_e_GeV = model->MakeField<double>("mu_e_GeV");
    reco.shower_e_GeV = model->MakeField<double>("shower_e_GeV");
    reco.sum_hadtrack_e_GeV = model->MakeField<double>("sum_hadtrack_e_GeV");

    truth.category = model->MakeField<int>("true_category");
    truth.cctopology = model->MakeField<int>("true_cctopology");
    truth.is_cc = model->MakeField<bool>("true_is_cc");
    truth.mode = model->MakeField<int>("true_mode");
    truth.vtx_pos_cm =
        model->MakeField<ROOT::Math::XYZPoint>("true_vtx_pos_cm");
    truth.pdgnu = model->MakeField<int>("true_pdgnu");
    truth.pdgnu_original = model->MakeField<int>("true_pdgnu_original");
    truth.enu_GeV = model->MakeField<double>("true_enu_GeV");
    truth.param_responses = model->MakeField<std::vector<std::vector<float>>>(
        "true_param_responses");
  }
};

} // namespace mach3winter26
