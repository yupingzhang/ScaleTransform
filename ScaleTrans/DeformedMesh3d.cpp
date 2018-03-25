#include "DeformedMesh3d.h"
#include <iostream>

using namespace std;

DeformedMesh3d::DeformedMesh3d(int numTri)
{
    // cout << "create Deformed Mesh3d..." << endl;
    triNum = numTri;
    
    _U.resize(triNum*3, 3);
    _V.resize(triNum*3, 3); 
    _s.resize(triNum, 3);

}

// input scale: (2.0, 1.0, 0.5)
void DeformedMesh3d::initDeformedState(float s1, float s2, float s3)
{
    for (int t = 0; t < triNum; ++t)
    {
        _s.row(t) = Eigen::Vector3d(s1, s2, s3);      

        _U.block<3,3>(3*t, 0) = Eigen::Matrix3d::Identity();
        _V.block<3,3>(3*t, 0) = Eigen::Matrix3d::Identity();

        theta_u.push_back(0);
        theta_v.push_back(0);
    }
}