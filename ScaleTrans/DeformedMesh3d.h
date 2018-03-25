// DeformedMesh3d

#ifndef SCALETRANS_DEFORMEDMESH3D_H
#define SCALETRANS_DEFORMEDMESH3D_H


#include <Eigen/Dense>
#include <vector>

class DeformedMesh3d
{
public:
	int triNum;
	Eigen::MatrixXd _U, _V, _s;   // left scaling, right scaling, singular values
    std::vector<double> theta_u, theta_v;  // rotation angle

	DeformedMesh3d(int numTris);
	~DeformedMesh3d() {};

	void initDeformedState(float s1, float s2, float s3);
};

#endif