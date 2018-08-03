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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  //Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  //Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  n_x_ = 5;

  n_aug_ = 7;

  n_z_radar_ = 3;

  n_z_lidar_ = 2;

  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

  // initial state vector
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);

  R_radar_ = MatrixXd::Zero(n_z_radar_,n_z_radar_);
  R_radar_(0,0) = std_radr_ * std_radr_;
  R_radar_(1,1) = std_radphi_ * std_radphi_;
  R_radar_(2,2) = std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd::Zero(n_z_lidar_,n_z_lidar_);
  R_lidar_(0,0) = std_laspx_ * std_laspx_;
  R_lidar_(1,1) = std_laspy_ * std_laspy_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if(!is_initialized_){
    cout<<"UKF:" <<endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      double py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      x_ << px, py, 0, 0, 0;
        /*double rho = meas_package.raw_measurements_(0);
        double phi = meas_package.raw_measurements_(1);
        double rhodot = meas_package.raw_measurements_(2);*/

        // polar to cartesian - r * cos(angle) for x and r * sin(angle) for y
        // ***** Middle value for 'v' can be tuned *****
        //x_ << rho * cos(phi), rho * sin(phi), 4, rhodot * cos(phi), rhodot * sin(phi);
        /*P_ << std_radr_*std_radr_, 0, 0, 0, 0,
                0, std_radr_*std_radr_, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, std_radphi_, 0,
                0, 0, 0, 0, std_radphi_;*/
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        //x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 4, 0.5, 0.0;

        //state covariance matrix
        //***** values can be tuned *****
        /*P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
                0, std_laspy_*std_laspy_, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;*/
    }

    P_<<0.15, 0, 0, 0, 0,
        0, 0.15, 0, 0, 0,
        0, 0, 0.15, 0, 0,
        0, 0, 0, 0.15, 0,
        0, 0, 0, 0, 0.15;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    cout<<"UKF: init done" <<endl;
    return;
  }

  double delta = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      //set measurement dimension, radar can measure r, phi, and r_dot
      UpdateRadar(meas_package);
  } else {
      //lidar can measure px, py
      UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    int n = 2 * n_aug_ + 1;
    lambda_ = 3 - n_aug_;

    weights_ = VectorXd::Zero(n);
    double weight = lambda_/(lambda_+ n_aug_);
    weights_(0) = weight;
    for (int i = 1; i < n; i++) {
        weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
    //create augmented mean vector
    VectorXd x_aug_ = VectorXd::Zero(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);

    //create sigma point matrix
    MatrixXd Xsig_aug_ = MatrixXd::Zero(n_aug_, n);


    //create augmented mean state
    x_aug_.head(n_x_) = x_;

    //create augmented covariance matrix
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;
    P_aug_(n_aug_ - 1, n_aug_ - 1) = std_yawdd_ * std_yawdd_;
    P_aug_(n_aug_ -2, n_aug_ -2) = std_a_ *std_a_;

    //create square root matrix
    MatrixXd P_aug_sqrt = P_aug_.llt().matrixL();

    //create augmented sigma points
    MatrixXd part1 = sqrt(n_aug_ + lambda_) * P_aug_sqrt;

    Xsig_aug_.col(0) = x_aug_;
    for(int i = 0; i < n_aug_; i++){
      Xsig_aug_.col(i+1) = x_aug_ + part1.col(i);
      Xsig_aug_.col(i + n_aug_ + 1) = x_aug_ - part1.col(i);
    }

    VectorXd c =  VectorXd(n_x_);
    VectorXd det =  VectorXd(n_x_);

    for(int i = 0; i < n; i++){
          //predict sigma points
          double v = Xsig_aug_(2,i);
          double psi = Xsig_aug_(3,i);
          double psi_dot = Xsig_aug_(4,i);
          double noise_a =  Xsig_aug_(5,i);
          double noise_psi =  Xsig_aug_(6,i);

          //deterministic part when psi dot not 0
          if(psi_dot != 0){
            det << (v/psi_dot) * (sin(psi+ (psi_dot * delta_t)) - sin(psi)),
                    (v/psi_dot) * (-cos(psi + (psi_dot * delta_t)) + cos(psi)),
                    0,
                    delta_t * psi_dot,
                    0;
          }else{
            //deterministic part when psi dot = 0
            det << delta_t * v * cos(psi),
                    delta_t * v * sin(psi),
                    0,
                    0,
                    0;
          }

          //stochastic part
          c << 0.5 * delta_t * delta_t * cos(psi) * noise_a,
                  0.5 * delta_t * delta_t * sin(psi) * noise_a,
                  delta_t * noise_a,
                  0.5 * delta_t * delta_t * noise_psi,
                  delta_t * noise_psi;

          Xsig_pred_.col(i) = Xsig_aug_.col(i).head(n_x_) + det + c;
      }

    //create vector for predicted state
    VectorXd x_pred = VectorXd::Zero(n_x_);

    //create covariance matrix for prediction
    MatrixXd P_pred = MatrixXd::Zero(n_x_, n_x_);

    //predict state mean
    for(int i = 0; i < n; i++){
        x_pred += weights_(i) * Xsig_pred_.col(i);
    }
  //predict state covariance matrix
    for(int i = 0; i < n; i++){
        //predict state covariance matrix
        VectorXd x_diff = Xsig_pred_.col(i) - x_pred;

        //normalize angles
        if (x_diff(3) > M_PI) {
            x_diff(3) -= 2. * M_PI;
        } else if (x_diff(3) < -M_PI) {
            x_diff(3) += 2. * M_PI;
        }

        P_pred +=  (weights_(i) * x_diff * x_diff.transpose());
    }

    x_ = x_pred;
    P_ = P_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  n_z_ = n_z_lidar_;
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd::Zero(n_z_, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z_);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z_,n_z_);

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

    MatrixXd R = R_lidar_;

    VectorXd z_ = VectorXd::Zero(n_z_);

    z_ = meas_package.raw_measurements_;

    //transform sigma points into measurement space
    for(int i =0; i <  2 * n_aug_ + 1; i++){
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        Zsig_.col(i) << px, py;
    }

    //calculate mean predicted measurement
    for(int i =0; i <  2 * n_aug_ + 1; i++){
        z_pred += (weights_(i) * Zsig_.col(i));
    }


    //calculate innovation covariance matrix S
    for(int i = 0; i < 2 * n_aug_ + 1; i++){
        VectorXd diff = (Zsig_.col(i) - z_pred);
        S +=(weights_(i) * diff * diff.transpose());
    }

    S += R;

    //calculate cross correlation matrix
    for(int i =0; i < 2 * n_aug_ + 1; i++){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        VectorXd Zdiff = Zsig_.col(i) - z_pred;
        //normalize angles
        if (x_diff(3) > M_PI) {
            x_diff(3) -= 2. * M_PI;
        } else if (x_diff(3) < -M_PI) {
            x_diff(3) += 2. * M_PI;
        }

        if(Zdiff(1)< - M_PI ){
            Zdiff(1) += 2 * M_PI;
        }else if (Zdiff(1) >  M_PI ){
            Zdiff(1) -= 2 * M_PI;
        }
        Tc += (weights_(i)* x_diff *  Zdiff.transpose());
    }

    //calculate Kalman gain K;
    MatrixXd k_gain = Tc * S.inverse();
    //update state mean and covariance matrix
    double temp = x_(0);

    VectorXd z_diff = z_ - z_pred;

    //normalize angles
    if (z_diff(1) > M_PI) {
        z_diff(1) -= 2. * M_PI;
    } else if (z_diff(1) < -M_PI) {
        z_diff(1) += 2. * M_PI;
    }
    x_ += (k_gain *z_diff);

    P_ -= k_gain * S * k_gain.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    n_z_ = n_z_radar_;
    //create matrix for sigma points in measurement space
    Zsig_ = MatrixXd::Zero(n_z_, 2 * n_aug_ + 1);

    //mean predicted measurement
    VectorXd z_pred = VectorXd::Zero(n_z_);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd::Zero(n_z_,n_z_);

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

    MatrixXd R = R_radar_;

    VectorXd z_ = VectorXd::Zero(n_z_);

    z_ = meas_package.raw_measurements_;

    //transform sigma points into measurement space
    for(int i =0; i <  2 * n_aug_ + 1; i++){
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double psi = Xsig_pred_(3,i);


        Zsig_.col(i) << sqrt((px * px) + (py *py)),
                atan2(py,px),
                ((px * cos(psi) * v) + (py * sin(psi) * v))/sqrt(px * px + py *py);

    }

    //calculate mean predicted measurement
    for(int i =0; i <  2 * n_aug_ + 1; i++){
        z_pred += (weights_(i) * Zsig_.col(i));
    }


    //calculate innovation covariance matrix S
    for(int i = 0; i < 2 * n_aug_ + 1; i++){
        VectorXd diff = (Zsig_.col(i) - z_pred);
        //normalize phi
        if(diff(1)< - M_PI ){
            diff(1) += 2 * M_PI;
        }else if (diff(1) >  M_PI ){
            diff(1) -= 2 * M_PI;
        }

        S +=  (weights_(i) * diff * diff.transpose());
    }

    S +=  R;

    //calculate cross correlation matrix
    for(int i =0; i < 2 * n_aug_ + 1; i++){
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        VectorXd Zdiff = Zsig_.col(i) - z_pred;
        //normalize angles
        if (x_diff(3) > M_PI) {
            x_diff(3) -= 2. * M_PI;
        } else if (x_diff(3) < -M_PI) {
            x_diff(3) += 2. * M_PI;
        }
        if(Zdiff(1)< - M_PI ){
            Zdiff(1) += 2 * M_PI;
        }else if (Zdiff(1) >  M_PI ){
            Zdiff(1) -= 2 * M_PI;
        }

        Tc += (weights_(i)* x_diff *  Zdiff.transpose());

    }


    //calculate Kalman gain K;
    MatrixXd k_gain = Tc * S.inverse();
    //update state mean and covariance matrix
    double temp = x_(0);

    VectorXd z_diff = z_ - z_pred;

    //normalize angles
    if (z_diff(1) > M_PI) {
        z_diff(1) -= 2. * M_PI;
    } else if (z_diff(1) < -M_PI) {
        z_diff(1) += 2. * M_PI;
    }
    x_ += (k_gain *z_diff);

    P_ -=  k_gain * S * k_gain.transpose();

}

