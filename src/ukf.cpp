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

  // if true data will beinitialized
  is_initialized_ = false;

  // initial value of time
  time_us_ = 0;

  ///* state vector dimension
  n_x_ = 5;

  ///* augmented state vector dimension
  n_aug_ = 7;

  // lamda factor
  lambda_ = 3 - n_aug_;

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

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0,    0,     0,
        0, 1, 0,    0,     0,
        0, 0, 1000, 0,     0,
        0, 0, 0,    1000 , 0, 
        0, 0, 0,    0,     1000;

  ///* augmented state vector
  x_aug_ = VectorXd(n_aug_);

  ///* augmented state covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  ///* augmented sigma points matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  ///* Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0,                       std_laspy_ * std_laspy_;

  ///* measurement radar noise matrix
  R_radar_ = MatrixXd(3, 3);
  R_lidar_ << std_radr_ * std_radr_, 0,                         0,
              0,                     std_radphi_ * std_radphi_, 0,
              0,                     0,                         std_radrd_ * std_radrd_;


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

      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];

      x_ << px, py, 0, 0, 0;

    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      double rho     = meas_package.raw_measurements_[0];
      double phi     = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      x_ << rho * cos( phi ), rho * sin( phi ), 0, 0, 0;

    }
  
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    
    for (int i=0; i<2*n_aug_+1;i++) {

      if (i == 0)
        weights_(i) = lambda_ / (lambda_ + 2*n_aug_+1);
      else
        weights_(i) = 1 / ( 2 * (lambda_ + 2*n_aug_+1) );
    
    }

    /**
      return from initialization.
      No need to predict or update.
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
  Augmentation();

  /**
    Create sigma points  
  */
  CreateSigmaPoints();

  /**
    Predict sigma points
  */
  Prediction(dt_);

  /**
    Predict mean and covariance
  */
  PredictMeanAndCovariance();


















  
  
  
  
  

 

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /*double px;
  double py;
  double yaw;
  double vx;
  double vy;
  double rho;
  double phi;
  double rho_d;

  Z_mea_lidar_sig = MatrixXd(4, 2*n_aug_+1);
  Z_mea_radar_sig = MatrixXd(3, 2*n_aug_+1);

  for(int i=0; i<2*n_aug_+1; i++) {

    px  = Xsig_pred_(0, i);
    py  = Xsig_pred_(1, i);
    yaw = Xsig_pred_(3, i);
    vx  = Xsig_pred_(2, i) * cos( yaw );
    vy  = Xsig_pred_(2, i) * sin( yaw );

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      Z_mea_lidar_sig(0, i) = px;      
      Z_mea_lidar_sig(1, i) = py;
      Z_mea_lidar_sig(2, i) = vx;
      Z_mea_lidar_sig(3, i) = vy;

    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

      rho = sqrt( ( px * px ) + ( py * py ) );
      phi = atan2( py, px );
      rho_d = ( ( px * vx ) + ( py * vy) ) / rho;

      Z_radar_sig_pred(0, i) = rho;      
      Z_radar_sig_pred(1, i) = phi;
      Z_radar_sig_pred(2, i) = rho_d;

    }

  }


  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    
    z_pred = VectorXd(3);
    z_pred.fill(0.0);
    S = MatriXd(3, 3);
    S.fill(0.0);
    R_radar = MatrixXd(3, 3);
    R_radar << std_radr_ * std_radr_, 0,                          0,
               0,                     std_radphi_ * std_radphi_ , 0,
               0,                     0,                          std_radrd_ * std_radrd_;

    for (int i=0; i<2*n_aug_+1; i++) {


      z_pred = z_pred + weights(i) * Z_radar_sig_pred.col(i);
      VectorXd z_diff = Z_radar_sig_pred.col(i) - z_pred;

      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      S = S + weights(i) * z_diff * z_diff.transpose();

    }

    S = S + R_radar;

    UpdateLidar(meas_package);

  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    z_pred = VectorXd(2);
    z_pred.fill(0.0);
    S = MatriXd(2, 2);
    S.fill(0.0);
    R_lidar = MatrixXd(2, 2);
    R_lidar << std_laspx_ * std_laspx_, 0,                       
               0,                     std_laspy_ * std_laspy_ ;

    for (int i=0; i<2*n_aug_+1; i++) {

      z_pred = z_pred + weights(i) * Z_mea_radar_sig.col(i);
      VectorXd z_diff = Z_mea_radar_sig.col(i) - z_pred;

      S = S + weights(i) * z_diff * z_diff.transpose();

    }

    S = S + R_lidar;

  }

  UpdateRadar(meas_package);*/

}

void UKF::Augmentation() {

  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;

}

void UKF::CreateSigmaPoints() {
  
  Xsig_aug_.col(0) = x_aug_;

  MatrixXd L;

  for (int i=0; i<n_aug_; i++) {

    L = P_aug_.llt().matrixL();

    Xsig_aug_.col(i+1)        = x_aug_ + sqrt( lambda_ + n_aug_ ) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt( lambda_ + n_aug_ ) * L.col(i);

  }

}

/**
 * Predicts sigma points.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  
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
  Predict the state and the state covariance matrix
*/
void UKF::PredictMeanAndCovariance() {
  
  x_.fill(0.0);
  P_.fill(0.0);

  for (int i=0; i<2*n_aug_+1; i++) {

    /**
      Predicted mean state vector
    */
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);

    VectorXd x_d = Xsig_pred_.col(i) - x_;
    while(x_d(3) > M_PI)
      x_d(3) -= 2 * M_PI;
    while(x_d(3) < -M_PI)
      x_d(3) += 2 * M_PI;

    /**
      Predicted state covariance matrix
    */
    P_ = P_ + weights_(i) * x_d * x_d.transpose();

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

  /**
    cross correlation matrix
  */
  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0.0);

  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }

  for (int i=0; i<2*n_aug_+1; i++) {

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Z_mea_radar_sig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    while (x_diff(3)> M_PI) x_diff(1)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(1)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();

  }

  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

}
