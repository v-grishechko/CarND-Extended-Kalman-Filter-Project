#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

    /**
    * A helper method to calculate Jacobians.
    */
    MatrixXd CalculateJacobian(const VectorXd &x_state);


    /**
     * A helper method to convert polar coordinates to cartesian
     */

    VectorXd ToCartesian(const VectorXd &polar_coord);

    /**
    * A helper method to convert cartesian coordinates to polar
    */
    VectorXd ToPolar(const VectorXd &cartesian_coord);

    float normalizeAngle(float angle);

};

#endif /* TOOLS_H_ */
