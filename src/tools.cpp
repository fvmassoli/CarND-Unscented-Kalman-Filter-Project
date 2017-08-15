#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  /**
  
  instanciate and initialize RMSE vector
  
  */
  VectorXd rmse(4);
  rmse.fill(0.0);

  /**
  
  the estimation vector size should not be 0
  the estimation vector size has to be same as the ground_truth size

  */
  if(estimations.size() != ground_truth.size() || estimations.size() == 0) {
    cout << "Invalid input. Cannot evaluate RMSE!" << endl;
    return rmse;
  }

  /**
  
  loop over the measurements

  */
  VectorXd tmp;
  for(int i=0; i<estimations.size(); i++) {
    tmp = estimations[i] - ground_truth[i];
    tmp = tmp.array()*tmp.array();
    rmse += tmp;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;

}