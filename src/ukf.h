#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  /**
   * Constructor
   */
  UKF(bool use_laser,
      bool use_radar,
      int n_x,
      int n_aug,
      double std_a,
      double std_yawdd,
      double std_laspx,
      double std_laspy,
      double std_radr,
      double std_radphi,
      double std_radrd);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(const MeasurementPackage &meas_package);

  inline VectorXd x() const { return x_; }
  inline double NIS_laser() const { return NIS_laser_; }
  inline double NIS_radar() const { return NIS_radar_; }
  
private:

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);
  MatrixXd GenerateAugmentedSigmaPts() const;
  void PredictSigmaPoints(MatrixXd &Xsig_pred,
                          const MatrixXd &Xsig_aug,
                          double delta_t) const;
  void PredictMeanAndCovariance(VectorXd &x,
                                MatrixXd &P,
                                const MatrixXd &Xsig_pred) const;

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage &meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage &meas_package);
  
  double NormalizeAngle(double angle) const;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;
  
  ///* State dimension - CTRV consists of [ px, py, v, psi, psi_dot ]
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Process noise standard deviation longitudinal acceleration in m2/s^4
  double std_a2_;

  ///* Process noise standard deviation yaw acceleration in rad2/s^4
  double std_yawdd2_;

  ///* Laser measurement noise standard deviation position1 in m2
  double std_laspx2_;

  ///* Laser measurement noise standard deviation position2 in m2
  double std_laspy2_;

  ///* Radar measurement noise standard deviation radius in m2
  double std_radr2_;

  ///* Radar measurement noise standard deviation angle in rad2
  double std_radphi2_;

  ///* Radar measurement noise standard deviation radius change in m2/s2
  double std_radrd2_ ;

  ///* time when the state is true, in us
  long long time_us_;

  ///* precalculated sqrt (lambda * n_aug)
  double lambda_sqrt_;
  
  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

};

#endif /* UKF_H */
