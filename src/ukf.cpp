#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double EPSILON = 1e-6;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF(bool use_laser,
         bool use_radar,
         int n_x,
         int n_aug,
         double std_a,
         double std_yawdd,
         double std_laspx,
         double std_laspy,
         double std_radr,
         double std_radphi,
         double std_radrd):
  is_initialized_(false),
  use_laser_(use_laser),
  use_radar_(use_radar),
  n_x_(n_x),
  n_aug_(n_aug),
  lambda_(3 - n_aug_),
  std_a2_(std_a * std_a),
  std_yawdd2_(std_yawdd * std_yawdd),
  std_laspx2_(std_laspx * std_laspx),
  std_laspy2_(std_laspy * std_laspy),
  std_radr2_(std_radr * std_radr),
  std_radphi2_(std_radphi * std_radphi),
  std_radrd2_(std_radrd * std_radrd),
  lambda_sqrt_(sqrt(lambda_ + n_aug_)),
  time_us_(0),
  x_(VectorXd::Zero(5)),
  P_(MatrixXd::Identity(5, 5)),
  Xsig_pred_(MatrixXd::Zero(5, 15)),
  NIS_radar_(0),
  NIS_laser_(0)
{
  weights_ = VectorXd(n_aug_ * 2 + 1);
  double w1 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = w1;
  double w2 = 0.5 / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); ++i) weights_(i) = w2;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &meas_package) {
  if (!is_initialized_)
  {
    switch (meas_package.sensor_type_)
    {      
    case MeasurementPackage::SensorType::LASER:
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      break;
    case MeasurementPackage::SensorType::RADAR:
    {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      x_ << px, py, 0, 0, 0;
      break;
    }
    default:
      std::cerr << "Unexpected SensorType - skipping measurement processing\n";
      return;
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  if ( (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER && !use_laser_) ||
       (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR && !use_radar_) )
    return;
  
  if (meas_package.timestamp_ > time_us_) {
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;
    Prediction(dt);
  }

  switch (meas_package.sensor_type_)
  {
  case MeasurementPackage::SensorType::LASER:
    UpdateLidar(meas_package);
    break;
  case MeasurementPackage::SensorType::RADAR:
    UpdateRadar(meas_package);
    break;
  default:
    std::cerr << "Unexpected SensorType - skipping measurement processing\n";
    break;
  }
}

MatrixXd UKF::GenerateAugmentedSigmaPts() const
{
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug.tail(n_aug_ - n_x_) << 0, 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.block(0, 0, n_x_, n_x_) << P_;
  P_aug.block(5, 5, 2, 2) << std_a2_, 0,
                             0, std_yawdd2_;

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd Al = lambda_sqrt_ * A;
  Xsig_aug.col(0) << x_aug;
  for (int i = 0; i < n_aug_; ++i) 
  {
    Xsig_aug.col(i + 1) << x_aug + Al.col(i);
    Xsig_aug.col(i + 1 + n_aug_) << x_aug - Al.col(i);
  }  

  return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd &Xsig_pred,
                             const MatrixXd &Xsig_aug,
                             double delta_t) const
{
  double delta_t2 = delta_t * delta_t;
  for (int i = 0; i < Xsig_pred.cols(); ++i)
  {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double psi = Xsig_aug(3, i);
    double psi_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);   // longitudinal acceleration noise
    double nu_psi = Xsig_aug(6, i); // yaw acceleration noise

    double cos_psi = cos(psi);
    double sin_psi = sin(psi);              
    VectorXd k1 = VectorXd::Zero(5);
    if (fabs(psi_dot) < EPSILON)
    {
      k1[0] = v * cos_psi * delta_t;
      k1[1] = v * sin_psi * delta_t;
    }
    else
    {
      k1[0] = (v / psi_dot) * (sin(psi + (psi_dot * delta_t)) - sin_psi);
      k1[1] = (v / psi_dot) * (-cos(psi + (psi_dot * delta_t)) + cos_psi);
      k1[3] = psi_dot * delta_t;
    }
    VectorXd noise = VectorXd(5);
    noise << 0.5 * delta_t2 * cos_psi * nu_a,
             0.5 * delta_t2 * sin_psi * nu_a,
             delta_t * nu_a,
             0.5 * delta_t2 * nu_psi,
             delta_t * nu_psi;
    
    VectorXd xk = VectorXd(5);
    xk << Xsig_aug.block(0, i, 5, 1);
    
    Xsig_pred.col(i) = xk + k1 + noise;
  }
}

double UKF::NormalizeAngle(double angle) const
{
  while (angle > M_PI) angle -= 2 * M_PI;
  while (angle <-M_PI) angle += 2 * M_PI;
  return angle;
}

void UKF::PredictMeanAndCovariance(VectorXd &x,
                                   MatrixXd &P,
                                   const MatrixXd &Xsig_pred) const
{
  for (int i = 0; i < n_x_; ++i) {
    x(i) = (Xsig_pred.row(i) * weights_).sum();
  }

  P.fill(0.0);
  for (int i = 0; i < Xsig_pred.cols(); i++) {
    VectorXd x_diff = Xsig_pred.col(i) - x;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }
}
  
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug = GenerateAugmentedSigmaPts();
  PredictSigmaPoints(Xsig_pred_, Xsig_aug, delta_t);
  PredictMeanAndCovariance(x_, P_, Xsig_pred_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage &meas_package) {
  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);  
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    Zsig.col(i) << px, py;
  }
  
  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);   
  for (int i = 0; i < z_pred.rows(); ++i) {
    z_pred(i) = (Zsig.row(i) * weights_).sum();
  }

  //calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S << std_laspx2_, 0,
       0, std_laspy2_;
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose() ;
  }

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd::Zero(n_x_, z_pred.rows());
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose() ;
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  x_ = x_ + K * z_diff;
  x_[3] = NormalizeAngle(x_[3]);
  P_ = P_ - K * S * K.transpose();

  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage &meas_package) {
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z); 
  MatrixXd S = MatrixXd(n_z, n_z);

  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double psi = Xsig_pred_(3, i);
    double rho = sqrt(px*px + py*py);
    double phi = atan2(py, px);
    if (fabs(rho) < EPSILON)
    {
      std::cerr << "UpdateRadar method failed, rho = 0. Skipping...\n";
      return;
    }
    double rho_dot = ((px * cos(psi) * v) + (py * sin(psi) * v)) / rho;
    Zsig.col(i) << rho, phi, rho_dot;
  }
  
  //calculate mean predicted measurement
  for (int i = 0; i < z_pred.rows(); ++i) {
    z_pred(i) = (Zsig.row(i) * weights_).sum();
  }

  //calculate measurement covariance matrix S
  S << std_radr2_, 0, 0,
       0, std_radphi2_, 0,
       0, 0, std_radrd2_;
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose() ;
  }

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd::Zero(n_x_, z_pred.rows());
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose() ;
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  z_diff(1) = NormalizeAngle(z_diff(1));
  x_ = x_ + K * z_diff;
  x_[3] = NormalizeAngle(x_[3]);
  P_ = P_ - K * S * K.transpose();

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
