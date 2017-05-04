
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"
#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void import_measurements(const string &input_file_path,
                         vector<MeasurementPackage> &measurement_pack_list,
                         vector<GroundTruthPackage> &gt_pack_list)
{
  ifstream in_file(input_file_path.c_str(), ifstream::in);
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << input_file_path << endl;
    exit(EXIT_FAILURE);
  }

  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // laser measurement
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // radar measurement
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  in_file.close();
}

int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string input_file_path = argv[1];
  string out_file_path = argv[2];

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;
  import_measurements(input_file_path, measurement_pack_list, gt_pack_list);
  
  UKF ukf(true, true, 5, 7,
          0.37, 0.56,
          0.15, 0.15,
          0.3, 0.03, 0.3);

  ofstream out_file(out_file_path.c_str(), ofstream::out);
  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_file_path << endl;
    exit(EXIT_FAILURE);
  }
  
  out_file << "time_stamp" << "\t";  
  out_file << "px_state" << "\t";
  out_file << "py_state" << "\t";
  out_file << "v_state" << "\t";
  out_file << "yaw_angle_state" << "\t";
  out_file << "yaw_rate_state" << "\t";
  out_file << "sensor_type" << "\t";
  out_file << "NIS" << "\t";  
  out_file << "px_measured" << "\t";
  out_file << "py_measured" << "\t";
  out_file << "px_ground_truth" << "\t";
  out_file << "py_ground_truth" << "\t";
  out_file << "vx_ground_truth" << "\t";
  out_file << "vy_ground_truth" << "\n";

  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;
  vector<double> nis_laser;
  vector<double> nis_radar;
  for (size_t k = 0; k < measurement_pack_list.size(); ++k) {
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // timestamp
    out_file << measurement_pack_list[k].timestamp_ << "\t"; // pos1 - est

    // output the state vector
    out_file << ukf.x()(0) << "\t"; // pos1 - est
    out_file << ukf.x()(1) << "\t"; // pos2 - est
    out_file << ukf.x()(2) << "\t"; // vel_abs -est
    out_file << ukf.x()(3) << "\t"; // yaw_angle -est
    out_file << ukf.x()(4) << "\t"; // yaw_rate -est

    // output lidar and radar specific data
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      out_file << "lidar" << "\t";
      out_file << ukf.NIS_laser() << "\t";
      out_file << measurement_pack_list[k].raw_measurements_(0) << "\t";
      out_file << measurement_pack_list[k].raw_measurements_(1) << "\t";
      nis_laser.push_back(ukf.NIS_laser());
    }
    else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      out_file << "radar" << "\t";
      out_file << ukf.NIS_radar() << "\t";
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file << ro * cos(phi) << "\t"; // px measurement
      out_file << ro * sin(phi) << "\t"; // py measurement
      nis_radar.push_back(ukf.NIS_radar());
    }

    // output the ground truth
    out_file << gt_pack_list[k].gt_values_(0) << "\t";
    out_file << gt_pack_list[k].gt_values_(1) << "\t";
    out_file << gt_pack_list[k].gt_values_(2) << "\t";
    out_file << gt_pack_list[k].gt_values_(3) << "\n";

    // convert ukf x vector to cartesian to compare to ground truth
    VectorXd ukf_x_cartesian_ = VectorXd(4);

    float x_estimate_ = ukf.x()(0);
    float y_estimate_ = ukf.x()(1);
    float vx_estimate_ = ukf.x()(2) * cos(ukf.x()(3));
    float vy_estimate_ = ukf.x()(2) * sin(ukf.x()(3));
    
    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;
    
    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }
  out_file.close();

  // compute the accuracy (RMSE)
  cout << "RMSE" << endl << Tools::CalculateRMSE(estimations, ground_truth) << endl;

  // Output NIS based on Chi-squared distribution
  int laser_in_range = 0;
  for (auto n : nis_laser) if (n > 0.103 && n < 5.991) laser_in_range++;
  cout << "NIS laser in range (expecting ~0.9) = " << float(laser_in_range) / nis_laser.size() << endl;
  int radar_in_range = 0;
  for (auto n : nis_radar) if (n > 0.352 && n < 7.815) radar_in_range++;
  cout << "NIS radar in range (expecting ~0.9) = " << float(radar_in_range) / nis_radar.size() << endl;
  cout << "Done!" << endl;
  return 0;
}
