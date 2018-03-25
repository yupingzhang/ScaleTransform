#pragma once

#include <Eigen/Dense>
#include <Eigen/Core>
#include <igl/writeOBJ.h>
#include <string>
#include <sstream>
#include <iomanip>

#include "DeformableMesh.h"


class TriangleMesh : public DeformableMesh
{
public:
    TriangleMesh(Eigen::MatrixXd V, Eigen::MatrixXd V1, Eigen::MatrixXd V2, Eigen::MatrixXi F);
    TriangleMesh(Eigen::MatrixXd V, Eigen::MatrixXi F);
    ~TriangleMesh();
  
    void initMesh(); 
	void initMeshDeform2d();
	void initMeshDeform3d();

    void addDeformationState2d(float p, Eigen::MatrixXd* deformation);
    void addDeformationState3d(float p, Eigen::MatrixXd* deformation);
    void recoverMesh2d(int i, Eigen::MatrixXd* deformGrad, Eigen::MatrixXd* deformMesh);
    void recoverMesh3d(int i, Eigen::MatrixXd* deformGrad, Eigen::MatrixXd* deformMesh);

};

