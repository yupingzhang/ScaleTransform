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

    _s.resize(triNum, 2);

    R.resize(triNum*2, 2);
    T.resize(triNum*2, 2);

    theta_u.resize(triNum);
    theta_v.resize(triNum);

    F = deformgrad;
}

// compute SVD
void DeformationGradient2d::svd()
{
    cout << "compute SVD..." << endl;     // compute svd per component
    for (int t = 0; t < triNum; ++t)
    {
    	Eigen::JacobiSVD<Eigen::MatrixXd> svd(F.block<2,2>(2*t, 0), Eigen::ComputeFullU | Eigen::ComputeFullV);

    	_s.row(t) = svd.singularValues();      

    	_U.block<2,2>(2*t, 0) = svd.matrixU();
    	_V.block<2,2>(2*t, 0) = svd.matrixV();

        Eigen::Rotation2Dd r1;
        r1.fromRotationMatrix(_U.block<2,2>(2*t, 0));
        theta_u.push_back(r1.angle());

        Eigen::Rotation2Dd r2;
        r2.fromRotationMatrix(_V.block<2,2>(2*t, 0));
        theta_v.push_back(r2.angle());

    }
	
}
