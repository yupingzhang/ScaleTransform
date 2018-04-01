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
void DeformedMesh3d::initDeformedState(double s1, double s2, double s3, double u, double v)
{
    for (int t = 0; t < triNum; ++t)
    {
        _s.row(t) = Eigen::Vector3d(s1, s2, s3);      

        // quaternion: 0+2i-j-3k
        vector<float> xaxis{1, 0, 0};
        vector<float> yaxis{0, 1, 0};
        float cos_u = cos(0.5*u);
        float sin_u = sin(0.5*u);
        float cos_v = cos(0.5*v);
        float sin_v = sin(0.5*v);

        quaternion_u = Eigen::Quaterniond(cos_u, xaxis[0] * sin_u, xaxis[1] * sin_u, xaxis[2] * sin_u);
        quaternion_v = Eigen::Quaterniond(cos_v, yaxis[0] * sin_v, yaxis[1] * sin_v, yaxis[2] * sin_v);

        _U.block<3,3>(3*t, 0) = quaternion_u.toRotationMatrix();
        _V.block<3,3>(3*t, 0) = quaternion_v.toRotationMatrix();

    }
}
