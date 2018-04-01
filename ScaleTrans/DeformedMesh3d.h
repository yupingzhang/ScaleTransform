// DeformedMesh3d

#ifndef SCALETRANS_DEFORMEDMESH3D_H
#define SCALETRANS_DEFORMEDMESH3D_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>


class DeformedMesh3d
{
public:
	int triNum;
	Eigen::MatrixXd F, _U, _V, _s;   // left scaling, right scaling, singular values
    Eigen::Quaterniond quaternion_u, quaternion_v;  // rotation angle

	DeformedMesh3d(int numTris);
	~DeformedMesh3d() {};

	void initDeformedState(double s1, double s2, double s3, double u, double v);
};

#endif