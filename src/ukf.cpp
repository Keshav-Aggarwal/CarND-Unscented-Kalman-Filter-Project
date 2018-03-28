#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;
  lambda_ = 3+n_x_;
  n_sig_ = n_aug_*2 + 1;
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  weights_ = VectorXd(n_aug_*2+1);
  P_aug_  = MatrixXd(n_aug_, n_aug_);
  x_.fill(0);

  MatrixXd R_laser_ = MatrixXd(2,2);
  R_laser_.fill(0.0);
  R_laser_(0,0) = std_laspx_*std_laspx_;
  R_laser_(1,1) = std_laspy_*std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  if(!is_initialized_){
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_<<meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],0,0,0;
    }
    else{
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rhodot = meas_package.raw_measurements_[2];

      float px = cos(phi) * rho;
      float py = sin(phi) * rho;
      float vx = cos(phi) * rhodot;
      float vy = sin(phi) * rhodot;
      float v = sqrt(vx*vx + vy*vy);

      x_<<px, py, v, 0, 0;
    }
    // Initialize weights
    weights_(0) = lambda_/(lambda_ + n_aug_);
    for(int i = 1;i<n_sig_; i++){
      weights_(i) = 0.5 / (lambda_ + n_aug_);
    }

    // Initialize covariance
    P_.fill(0);
    P_(0,0) = P_(1,1) = P_(2,2) = P_(3,3) = P_(4,4) = 10;

    previousTimestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  float dt =  (meas_package.timestamp_ - previousTimestamp_)/1000000.0;
  previousTimestamp_ = meas_package.timestamp_;

  Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);
  }
  else{
    //UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  double delta_t2 = delta_t*delta_t;
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  //Fill augmented matrix
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;

  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_);

  for(int i =0; i<=n_aug_*2;i++)
  {
      double px = Xsig_aug(0,i);
      double py = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double yaw = Xsig_aug(3,i);
      double yawd = Xsig_aug(4,i);
      double nu_a = Xsig_aug(5,i);
      double nu_yawdd = Xsig_aug(6,i);
      //avoid division by zero
      double px_p = 0;
      double py_p = 0;
      if(fabs(yawd) > 0.0001){
         px_p = px + (v/yawd)*(sin(yaw+yawd*delta_t) - sin(yaw));
         py_p = py + (v/yawd)*(-cos(yaw+yawd*delta_t) + cos(yaw));
      }
      else{
         px_p = px + (v*cos(yaw)*delta_t); 
         py_p = py + (v*sin(yaw)*delta_t); 
      }
      double v_p = v;
      double yaw_p = yaw + yawd * delta_t;
      double yawd_p = yawd;
      
      px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
      py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
      v_p = v_p + nu_a * delta_t;
      
      yaw_p = yaw_p + 0.5 * delta_t*delta_t*nu_yawdd;
      yawd_p = yawd_p + delta_t * nu_yawdd;
      //write predicted sigma points into right column
      Xsig_pred(0,i) = px_p;
      Xsig_pred(1,i) = py_p;
      Xsig_pred(2,i) = v_p;
      Xsig_pred(3,i) = yaw_p;
      Xsig_pred(4,i) = yawd_p;
  }
  x_.fill(0.0);
   for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred.col(i);

  }
  
  //predict state covariance matrix
  P_.fill(0);

  for(int i=0;i<=n_aug_*2;i++){
      VectorXd x_diff = Xsig_pred.col(i) - x_;
      x_diff(3) = atan2(cos(x_diff(3)), sin(x_diff(3)));
      P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z = 2;
  MatrixXd zsig = Xsig_pred_.block(0, 0, n_z, n_sig_);
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  for (int i=0; i < n_sig_  ; i++) {
      z_pred = z_pred + weights_(i) * zsig.col(i);
  }

  MatrixXd S = MatrixXd(n_z, n_z);
  
  S.fill(0.0);
  for(int i =0;i<n_sig_;i++){
     VectorXd z_diff = zsig.col(i) - z_pred;
     S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0 ,
  0, std_laspy_*std_laspy_;

  
  S = S + R;
  
  MatrixXd T = MatrixXd(n_x_, n_z);
  T.fill(0.0);

  //* Cross correlation matrix
  for(int i=0; i<n_sig_; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    T = T + (weights_(i) * (x_diff)) * (zsig.col(i) - z_pred).transpose();
  }

  //Kalman Gain
  MatrixXd K = T * S.inverse();
  VectorXd z = meas_package.raw_measurements_;

  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {


}
