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
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd RMSE = VectorXd(4);
  RMSE<<0,0,0,0;
  for(int i=0;i<estimations.size();i++){
    VectorXd temp = estimations[i] - ground_truth[i];
    temp = temp.array() * temp.array();
    RMSE += temp;
  }
  //CALCULATING MEAN BY DIVING THE SQUARED DIFFERENCE BY NUMBER OF ELEMNETS
  RMSE = RMSE / estimations.size();
  //CALCULATING THE SQUARE ROOT OF THE ELEMENT
  RMSE = RMSE.array().sqrt();
  
  return RMSE;
}