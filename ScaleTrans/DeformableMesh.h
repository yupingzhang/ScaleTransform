//
// Created by Yuping Zhang on 3/2/18.
//

#ifndef SCALETRANS_DEFORMABLEMESH_H
#define SCALETRANS_DEFORMABLEMESH_H

#include <vector>
#include "DeformationGradient2d.h"

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

    Eigen::MatrixXd Bs;   // base matrix. (2, t)
    Eigen::MatrixXd X_st;   // deformation matrix
    Eigen::MatrixXd X_ed;   // deformation matrix

    // new vertex positions after interpolation
    vector< Eigen::MatrixXd> predictedVertices;


    DeformableMesh() {};
    virtual ~DeformableMesh() {};

};


#endif //SCALETRANS_DEFORMABLEMESH_H
