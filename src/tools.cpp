#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if ((estimations.empty()) || (estimations.size() != ground_truth.size())) {
    std::cout << "RMSE cannot be computed. " 
      "Check estimation and ground_truth vectors." << std::endl;
    return rmse;
  };

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations.at(i) - ground_truth.at(i);
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate the mean
  int num_samples = (int)estimations.size();
  rmse = rmse / num_samples;

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;  
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // compute help variables
  float px2 = px * px;
  float py2 = py * py;
  float rho2 = px2 + py2;
  float rho = sqrt(rho2);

  // check division by zero
  if (std::isinf(1/rho) || std::isinf(1/rho2)) {
    std::cout << "Division by zero. Hj cannot be computed." << std::endl;
  } else {
    // compute the Jacobian matrix
    float hj00 = px / rho;
    float hj01 = py / rho;
    float hj10 = -py / rho2;
    float hj11 = px / rho2;
    float hj20 = py * (vx * py - vy * px) / pow(rho2, 1.5);
    float hj21 = px * (vy * px - vx * py) / pow(rho2, 1.5);
    float hj22 = px / rho;
    float hj23 = py / rho;
    
    Hj << hj00, hj01, 0, 0,
          hj10, hj11, 0, 0,
          hj20, hj21, hj22, hj23;
  };
  return Hj;
}