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

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  time_us_ = 0;

  ///* Weights of sigma points
  VectorXd weights_;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

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
