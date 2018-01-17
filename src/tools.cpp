#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size()
        || estimations.size() == 0) {
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for (unsigned int i = 0; i < estimations.size(); ++i) {

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array() * residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
    MatrixXd Hj(3, 4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px * px + py * py;
    float c2 = sqrt(c1);
    float c3 = (c1 * c2);

    //check division by zero
    if (fabs(c1) < 0.0001) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }

    //compute the Jacobian matrix
    Hj << (px / c2), (py / c2), 0, 0,
            -(py / c1), (px / c1), 0, 0,
            py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;

    return Hj;
}

VectorXd Tools::ToCartesian(const VectorXd &polar_coord) {
    VectorXd cartesian(4);

    float px = polar_coord[0] * cos(polar_coord[1]);
    float py = polar_coord[1] * cos(polar_coord[2]);

    cartesian << px, py, 0, 0;

    return cartesian;
}

VectorXd Tools::ToPolar(const VectorXd &cartesian_coord) {
    VectorXd polar(3);

    float px = cartesian_coord[0];
    float py = cartesian_coord[1];
    float vx = cartesian_coord[2];
    float vy = cartesian_coord[3];

    float rho = sqrt(px * px + py * py);
    float theta = atan2(py, px);
    float rho_dot = ((px * vx) + (py * vy)) / (sqrt(px * px + py * py));

    polar << rho, theta, rho_dot;

    return polar;
}

float Tools::normalizeAngle(float angle) {
    if (angle > 2 * M_PI) {
        angle -= 2 * M_PI;
    } else if (angle < -2 * M_PI) {
        angle += 2 * M_PI;
    }
    return angle;
}
