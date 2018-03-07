//
// Created by Yuping Zhang on 3/2/18.
//

#ifndef SCALETRANS_DEFORMATIONGRADIENT2D_H
#define SCALETRANS_DEFORMATIONGRADIENT2D_H

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>

/**!
 *
 */
class DeformationGradient2d {
public:
    int triNum;
    Eigen::MatrixXd F;   // deformation gradient
    Eigen::MatrixXd R, T;   // rotation and stretch
    // SVD
    Eigen::MatrixXd _U, _V, _Vt, _S;   // left scaling, right scaling, singular values

    // create from known state
    DeformationGradient2d(int numTris, Eigen::MatrixXd deformgrad);

    ~DeformationGradient2d() {};

    // compute SVD
    void svd();

};


#endif //SCALETRANS_DEFORMATIONGRADIENT2D_H
