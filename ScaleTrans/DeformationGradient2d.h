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
    Eigen::MatrixXd F;   // deformation gradient
    Eigen::MatrixXd R, T;   // rotation and stretch
    // SVD
    Eigen::MatrixXd _U, _V, _Vt;
    Eigen::VectorXd _s;   // left scaling, right scaling, singular values  

    // new state without initial variable
    DeformationGradient2d() {};
    // create from known state
    DeformationGradient2d(Eigen::MatrixXd deformgrad);

    ~DeformationGradient2d() {};

    // compute SVD
    void svd();

    // update rotation and stretch matrix
    void updateRT( Eigen::MatrixXd U,  Eigen::MatrixXd V,  Eigen::VectorXd s );

};


#endif //SCALETRANS_DEFORMATIONGRADIENT2D_H
