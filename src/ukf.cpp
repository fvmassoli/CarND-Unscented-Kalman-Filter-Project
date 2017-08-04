#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  time_us_ = 0;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* state vector dimension
  n_x_ = 5;

  ///* augmented state vector dimension
  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);;

  ///* augmented state vector
  x_aug_ = VectorXd(n_aug_);

  ///* augmented state covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  ///* augmented sigma points matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_) {

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      px_ = meas_package.raw_measurements_[0];
      py_ = meas_package.raw_measurements_[1];

      x_ << px_, py_, 0, 0, 0;

    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      rho_ = meas_package.raw_measurements_[0];
      phi_ = meas_package.raw_measurements_[1];
      rho_dot_ = meas_package.raw_measurements_[2];

      x_ << rho_*cos(phi_), rho_*sin(phi_), 0, 0, 0;

    }
  
    P_ << 1, 0, 0,    0,     0,
          0, 1, 0,    0,     0,
          0, 0, 1000, 0,     0,
          0, 0, 0,    1000 , 0, 
          0, 0, 0,    0,     1000;

    time_us_ = meas_package.timestamp_;
    
    Augmentation(x_, P_, std_a_, std_yawdd_, n_x_, is_initialized_);
    
    is_initialized_ = true;

    /**
      return from initialization
    */
    return;
  }


  /*****************************************************************************
     Prediction
   ****************************************************************************/

  dt_ = ( meas_package.timestamp_ - time_us_ ) / 1000000.0; //  in seconds
  time_us_ = meas_package.timestamp_;

  /**
    Augmentation step
  */
  Augmentation(x_, P_, std_a_, std_yawdd_, n_x_, is_initialized_);
  
  /**
    Create sigma points  
  */
  CreateSigmaPoints(x_aug_, P_aug_, lambda_, n_aug_);
  
  /**
    Predict sigma points
  */
  Prediction(dt_);


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

  }



}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  double px, px_p;
  double py, py_p;
  double v, v_p;
  double yaw, yaw_p;
  double yaw_d, yaw_d_p;
  double nu_a;
  double nu_yawdd;
  double delta_t2 = delta_t * delta_t;

  for (int i=0; i<2*n_aug_+1; i++) {

    px       = Xsig_aug_(0, i);
    py       = Xsig_aug_(1, i);
    v        = Xsig_aug_(2, i);
    yaw      = Xsig_aug_(3, i);
    yaw_d    = Xsig_aug_(4, i);
    nu_a     = Xsig_aug_(5, i);
    nu_yawdd = Xsig_aug_(6, i);

    /**
      Avoid division by 0
    */
    if (fabs(yaw) > 0.001) {
      px_p = px + ( v / yaw_d) * ( sin( yaw + yaw_d*delta_t ) - sin( yaw ) ) + 0.5 * delta_t2 * cos(yaw) * nu_a;
      py_p = py + ( v / yaw_d) * ( -cos( yaw + yaw_d*delta_t ) + cos( yaw ) ) + 0.5 * delta_t2 * sin(yaw) * nu_a;
    } else {
      px_p = v * cos(yaw);
      py_p = v * sin(yaw);
    }

    v_p = nu_a * delta_t;
    
    if (fabs(yaw) > 0.001) {
      yaw_d_p = yaw_d * delta_t;
    } else {
      yaw_d_p = 0;
    } 

    yaw_d_p = nu_yawdd * delta_t;

    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, 1) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yaw_d_p;

  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

void UKF::Augmentation(VectorXd x, MatrixXd P, float std_a, float std_yawdd, int n_x, bool is_initialized) {

  if (!is_initialized) {
   
    x_aug_.head(n_x) = x;
    x_aug_(5) = 0;
    x_aug_(6) = 0;

    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_, n_x_) = P;
    P_aug_(5, 5) = std_a*std_a;
    P_aug_(6, 6) = std_yawdd*std_yawdd;

  }

} 

void UKF::CreateSigmaPoints(VectorXd x_aug, MatrixXd P_aug, double lambda, int n_aug) {
  
  Xsig_aug_.col(0) = x_aug;
  MatrixXd L;

  for (int i=0; i<n_aug_; i++) {

    L = P_aug.llt().matrixL();

    Xsig_aug_.col(i+1)        = x_aug + sqrt( lambda + n_aug ) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt( lambda + n_aug ) * L.col(i);

  }

}

void UKF::PredictMeanAndCovariance() {

}
