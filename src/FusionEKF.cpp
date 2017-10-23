#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
             0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
             0, 0.0009, 0,
             0, 0, 0.09;

    /**
    TODO:
      * Finish initializing the FusionEKF.
      * Set the process and measurement noises
    */
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ <<  1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ <<    1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ <<    1, 0, 0, 0,
            0, 1, 0, 0;

    ekf_.R_laser_ = R_laser_;
    ekf_.R_radar_ = R_radar_;

    noise_ax_ = 9.0;
    noise_ay_ = 9.0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF initialization" << endl;

        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        previous_timestamp_ = measurement_pack.timestamp_;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            double const x = measurement_pack.raw_measurements_(0) * std::cos(measurement_pack.raw_measurements_(1));
            double const y = measurement_pack.raw_measurements_(0) * std::sin(measurement_pack.raw_measurements_(1));

            ekf_.x_ << x, y, 1, 1;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 1, 1;
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     TODO:
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    double const dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    double const dt2 = dt * dt;
    double const dt3 = dt2 * dt;
    double const dt4 = dt3 * dt;

    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt4 / 4 * noise_ax_, 0, dt3 / 2 * noise_ax_, 0,
            0, dt4 / 4 * noise_ay_, 0, dt3 / 2 * noise_ay_,
            dt3 / 2 * noise_ax_, 0, dt2*noise_ax_, 0,
            0, dt3 / 2 * noise_ay_, 0, dt2*noise_ay_;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        ekf_.Hj_ = tools.CalculateJacobian( ekf_.x_ );
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "t = " << measurement_pack.timestamp_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
}
