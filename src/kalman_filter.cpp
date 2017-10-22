#include "kalman_filter.h"
#include "math.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  P_ = F_ * P_ * F_.transpose() + Q_;
  x_ = F_ * x_; 
}

void KalmanFilter::Update(const VectorXd &z) {
  MatrixXd y = z - H_ * x_;
  MatrixXd S = H_* P_ * H_.transpose() + R_laser_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I  - K * H_) * P_;    
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

	double const px = x_(0);
	double const py = x_(1);
	double const vx = x_(2);
	double const vy = x_(3);

	VectorXd z_p(3);
	double const range = std::abs(std::complex<double>(px, py));
	double const phi = std::atan2(py,  px);
	z_p << range, phi, (px*vx + py*vy) / range;

	VectorXd y = z - z_p;
	if (y(1) > M_PI )
    {
		y(1) -= 2 * M_PI ;
	}
    else if (y(1) < -M_PI )
	{
        y(1) += M_PI ;
    }

	MatrixXd Ht = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Ht + R_radar_;
	MatrixXd K = P_ * Ht * S.inverse();

	x_ = x_ + K * y;
	MatrixXd I = MatrixXd::Identity( x_.size(),  x_.size() );
	P_ = (I - K * Hj_) * P_;
}
