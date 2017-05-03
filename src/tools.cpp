#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse = VectorXd::Zero(4);
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    std::cerr << "ERROR: mismatched sizes for estimations and ground truth vectors."
              << std::endl;
    return rmse;
  }

  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd diff = estimations[i] - ground_truth[i];
    VectorXd diff2 = diff.array() * diff.array();
    rmse = diff2 + rmse;
  }
  
  rmse = rmse.array() / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;  
}
