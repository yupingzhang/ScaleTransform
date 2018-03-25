//
// Created by Yuping Zhang on 3/2/18.
//

#ifndef SCALETRANS_DEFORMABLEMESH_H
#define SCALETRANS_DEFORMABLEMESH_H


#include "DeformationGradient2d.h"
#include "DeformedMesh3d.h"

using namespace std;


class DeformableMesh {
public:
    // vertex positions of initial rest pose
    Eigen::MatrixXd initialVertices;     // (n, 3)
    // vertex indices of each triangle
    Eigen::MatrixXi triangles;

    // vertex positions of current mesh
    Eigen::MatrixXd deformedVertices_st;
    Eigen::MatrixXd deformedVertices_ed;

    // deformed states
    DeformationGradient2d* deformedState_st;
    DeformationGradient2d* deformedState_ed;

    DeformedMesh3d* deformedState3d_st;
    DeformedMesh3d* deformedState3d_ed;

    Eigen::MatrixXd Bs;   // base dx

    Eigen::MatrixXd F_st;   // deformation matrix
    Eigen::MatrixXd F_ed;   // deformation matrix


    DeformableMesh() {};
    virtual ~DeformableMesh() {};

};


#endif //SCALETRANS_DEFORMABLEMESH_H
