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

void KalmanFilter::RecalculateParameters(const MatrixXd &R, const MatrixXd &H, const VectorXd &y)
{
    MatrixXd P_Ht = P_ * H.transpose();
    MatrixXd S = H * P_Ht + R;
    MatrixXd K = P_Ht * S.inverse();

    x_ = x_ + K * y;
    MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
    P_ = (I  - K * H) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {

    VectorXd y = z - H_ * x_;
    RecalculateParameters(R_laser_, H_, y);
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

    RecalculateParameters(R_radar_, Hj_, y);
}
