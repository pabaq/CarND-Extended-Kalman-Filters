#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;  // will be set to Hj_ or H_laser_ before calling update step
  R_ = R_in;  // will be set to R_radar_ or R_laser_ before calling update step
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */

  unsigned int z_size = z.size();  // 2 (laser measures px and py)
  unsigned int x_size = x_.size();  // 4 (px, py, vx, vy)
  MatrixXd S(z_size, z_size);  // (2, 2)
  MatrixXd K(x_size, z_size);  // (4, 2)
  MatrixXd I = MatrixXd::Identity(x_size, x_size);  // (4, 4)
  VectorXd z_pred(z_size);  // (2,)
  VectorXd y(z_size);  // (2,)
  
  z_pred = H_ * x_;  // predicted state in measurement space
  y = z - z_pred;  // error between measurement and predicted state (innovation)
  S = H_ * P_ * H_.transpose() + R_;  // helper matrix
  K = P_ * H_.transpose() * S.inverse();  // Kalman Gain
  
  // new estimate
  x_ += K * y;
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  unsigned int z_size = z.size();  // 3 (laser measures rho, phi and rhodot)
  unsigned int x_size = x_.size();  // 4 (px, py, vx, vy)
  MatrixXd S(z_size, z_size);  // (3, 3)
  MatrixXd K(x_size, z_size);  // (4, 3)
  MatrixXd I = MatrixXd::Identity(x_size, x_size);  // (4, 4)
  VectorXd z_pred(z_size);  // (3,)
  VectorXd y(z_size);  // (3,)

  // recover state parameters in cartesian coordinates
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  // convert state from cartesian into polar coordinates (h function)
  double rho = sqrt(px * px + py * py);
  double phi = atan2(py, px);  // in interval [-pi, pi]
  double rhodot = (px * vx + py * vy) / std::max(rho, 1e-5);  // ommit division by 0

  z_pred << rho, phi, rhodot;
  y = z - z_pred;  // error between measurement and predicted state (innovation)
  S = H_ * P_ * H_.transpose() + R_;  // helper matrix
  K = P_ * H_.transpose() * S.inverse();  // Kalman Gain

  // Adjust phi in y to lie in interval [-pi, pi]
  while (y(1) < M_PI) { y(1) += M_PI*2; }
  while (y(1) > M_PI) { y(1) -= M_PI*2; }

  // new estimate
  x_ += K * y;
  P_ = (I - K * H_) * P_;
}
