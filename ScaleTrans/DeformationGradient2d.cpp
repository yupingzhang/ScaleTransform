//
// Created by Yuping Zhang on 3/2/18.
//

#include "DeformationGradient2d.h"
#include <iostream>

using namespace std;

DeformationGradient2d::DeformationGradient2d(Eigen::MatrixXd deformgrad)
{
    cout << "create deformation gradient..." << endl;
    F = deformgrad;
}

// compute SVD
void DeformationGradient2d::svd()
{
    cout << "compute SVD..." << endl;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // Eigen::BDCSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);

    _s = svd.singularValues();
    _U = svd.matrixU();
    _V = svd.matrixV();
    _Vt = _V.transpose().conjugate();

    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << _U.rows() << " x " << _U.cols() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << _V.rows() << " x " << _V.cols() << endl;

    //debug
    cout << "_U.log() ===================================================================================== " << endl;
    cout << _U.log() << endl;
}


// update rotation and stretch matrix
void DeformationGradient2d::updateRT( Eigen::MatrixXd U,  Eigen::MatrixXd V,  Eigen::VectorXd s )
{
    cout << "update R & T..." << endl;
    // convert s (1d vector for singular values) to square matrix
    Eigen::MatrixXd S = s.asDiagonal();
    Eigen::MatrixXd Vt = V.transpose().conjugate();
    R = U * Vt;
    T = V * S * Vt;
}