#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    int size = estimations.size();
    VectorXd sum(4);
    VectorXd rmse(4);
    rmse << 0 ,0 ,0 ,0;
    sum << 0, 0, 0, 0;

    if(size != ground_truth.size() || size == 0){
        cerr << "Error: Estimation and ground truth vectors sizes are not equal or zero.";
    } else{
        for (int i = 0; i < size; ++i) {
            VectorXd residual = estimations[i] - ground_truth[i];
            residual = residual.array().pow(2);
            sum += residual;
        }
        sum /= size;
        rmse << sum.array().sqrt();
    }
    return rmse;

}