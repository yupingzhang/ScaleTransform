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
//    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
//    Eigen::BDCSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);

    _s = svd.singularValues();
    _U = svd.matrixU();
    _V = svd.matrixV();
    _Vt = _V.transpose().conjugate();

    cout << "Its singular values are:" << endl << svd.singularValues().rows() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << _U.rows() << " x " << _U.cols() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << _V.rows() << " x " << _V.cols() << endl;

}
