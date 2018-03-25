//
// Created by Yuping Zhang on 3/2/18.
//

#include "DeformationGradient2d.h"
#include <iostream>

using namespace std;

DeformationGradient2d::DeformationGradient2d(int numTri)
{
    cout << "create deformation gradient..." << endl;
    triNum = numTri;
    
    _U.resize(triNum*2, 2);
    _V.resize(triNum*2, 2); 

    _s.resize(triNum, 2);

    R.resize(triNum*2, 2);
    T.resize(triNum*2, 2);

    
}

// compute SVD
void DeformationGradient2d::svd(Eigen::MatrixXd deformgrad)
{
    F = deformgrad;

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

        // if (t == 10)
        // {
        //     cout << "_U: " << endl;
        //     cout << _U.block<2,2>(2*t, 0) << endl;
        //     cout << "_V: " << endl;
        //     cout << _V.block<2,2>(2*t, 0) << endl;
        //     cout << "_s:" << endl;
        //     cout << _s.row(t) << endl;

        //     cout << "angles: " << r1.angle() << " " << r2.angle() << endl;     // problem: 0 ~ 1.57
        // }
    }
}


/**
 * @brief      create example deformation state
 */
void DeformationGradient2d::initDeformedState(float scale1, float scale2)
{
    for (int t = 0; t < triNum; ++t)
    {
        _s.row(t) = Eigen::Vector2d(scale1, scale2);      

        _U.block<2,2>(2*t, 0) = Eigen::Matrix2d::Identity();
        _V.block<2,2>(2*t, 0) = Eigen::Matrix2d::Identity();

        Eigen::Rotation2Dd r1(0);
        theta_u.push_back(r1.angle());

        Eigen::Rotation2Dd r2(0);
        theta_v.push_back(r2.angle());

    }
}



