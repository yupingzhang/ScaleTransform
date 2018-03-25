//
// Created by Yuping Zhang on 3/2/18.
//

#ifndef SCALETRANS_DEFORMATIONGRADIENT2D_H
#define SCALETRANS_DEFORMATIONGRADIENT2D_H

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Geometry> 
// #include <unsupported/Eigen/MatrixFunctions>

#include <vector>

/**!
 *
 */
class DeformationGradient2d {
public:
    int triNum;
    Eigen::MatrixXd F;   // deformation gradient
    Eigen::MatrixXd R, T;   // rotation and stretch
    // SVD
    Eigen::MatrixXd _U, _V, _s;   // left scaling, right scaling, singular values
    std::vector<double> theta_u, theta_v;  // rotation angle

    // create from known state
    DeformationGradient2d(int numTris);

    ~DeformationGradient2d() {};

    // compute SVD
    void svd(Eigen::MatrixXd deformgrad);

    void initDeformedState(float scale1, float scale2);

};


#endif //SCALETRANS_DEFORMATIONGRADIENT2D_H
