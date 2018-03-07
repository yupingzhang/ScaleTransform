//
// Created by Yuping Zhang on 3/2/18.
//

#include "DeformationGradient2d.h"
#include <iostream>

using namespace std;

DeformationGradient2d::DeformationGradient2d(int numTri, Eigen::MatrixXd deformgrad)
{
    cout << "create deformation gradient..." << endl;
    triNum = numTri;
    
    F.resize(triNum*2, 2);
    _U.resize(triNum*2, 2);
    _V.resize(triNum*2, 2); 
    _Vt.resize(triNum*2, 2); 
    _S.resize(triNum*2, 2);
    R.resize(triNum*2, 2);
    T.resize(triNum*2, 2);

    F = deformgrad;
}

// compute SVD
void DeformationGradient2d::svd()
{
    cout << "compute SVD..." << endl;     // compute svd per component
    for (int t = 0; t < triNum; ++t)
    {
    	Eigen::JacobiSVD<Eigen::MatrixXd> svd(F.block<2,2>(2*t, 0), Eigen::ComputeFullU | Eigen::ComputeFullV);

    	_S.block<2,2>(2*t, 0) = svd.singularValues().asDiagonal();
    	_U.block<2,2>(2*t, 0) = svd.matrixU();
    	_V.block<2,2>(2*t, 0) = svd.matrixV();
    	_Vt.block<2,2>(2*t, 0) = _V.transpose().conjugate();
    }
	/*
	    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

	    _S = svd.singularValues().asDiagonal();
	    _U = svd.matrixU();
	    _V = svd.matrixV();
	    _Vt = _V.transpose().conjugate();

	    R = U * Vt;   
	    T = V * S * Vt;
	*/

//    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
//    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << _U.rows() << " x " << _U.cols() << endl;
//    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << _V.rows() << " x " << _V.cols() << endl;

}
