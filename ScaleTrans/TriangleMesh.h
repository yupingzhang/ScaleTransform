#pragma once

#include <Eigen/Dense>
#include <Eigen/Core>

#include "DeformableMesh.h"
#include "DeformationGradient2d.h"


class TriangleMesh : public DeformableMesh
{
public:
    TriangleMesh(Eigen::MatrixXd V, Eigen::MatrixXd V1, Eigen::MatrixXd V2, Eigen::MatrixXi F);
    ~TriangleMesh();
  
    void initMesh(); 
    void addDeformationState(float t);

};

