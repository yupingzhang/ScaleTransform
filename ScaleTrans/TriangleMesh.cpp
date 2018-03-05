#include "TriangleMesh.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace std;

// initialize rest vertices, current vertices, and triangle indices
TriangleMesh::TriangleMesh(Eigen::MatrixXd V, Eigen::MatrixXd V1, Eigen::MatrixXd V2, Eigen::MatrixXi F)
{
  initialVertices = Eigen::MatrixXd(V);
  deformedVertices_st = Eigen::MatrixXd(V1);
  deformedVertices_ed = Eigen::MatrixXd(V2);
  triangles = Eigen::MatrixXi(F);
}

TriangleMesh::~TriangleMesh()
{
    predictedVertices.clear();
}

/**
 * { init the base matrix and two deformation gradient }
 */
void TriangleMesh::initMesh()
{
    cout << "====== initialize the mesh data =====" << endl;
    // Compute the base matrix for deformation gradients
    int numTris = triangles.rows();    // get triangle number
    Bs.resize(2, 2*numTris);
    X_st.resize(2, 2*numTris);
    X_ed.resize(2, 2*numTris);

    for(int t=0;t<numTris;t++)
    {
        Eigen::Vector2d A,B,C;
        Eigen::Matrix2d V;

        A = initialVertices.row(triangles(t,0)).block<1,2>(0,0);
        B = initialVertices.row(triangles(t,1)).block<1,2>(0,0);
        C = initialVertices.row(triangles(t,2)).block<1,2>(0,0);

        V << A-C,B-C;
        Bs.block<2,2>(0,2*t) = V.inverse().cast<double>();

        ////// deformation start
        A = deformedVertices_st.row(triangles(t,0)).block<1,2>(0,0);
        B = deformedVertices_st.row(triangles(t,1)).block<1,2>(0,0);
        C = deformedVertices_st.row(triangles(t,2)).block<1,2>(0,0);
        
        V << A-C,B-C;
        X_st.block<2,2>(0,2*t) = V.cast<double>();

        ////// deformation end
        A = deformedVertices_ed.row(triangles(t,0)).block<1,2>(0,0);
        B = deformedVertices_ed.row(triangles(t,1)).block<1,2>(0,0);
        C = deformedVertices_ed.row(triangles(t,2)).block<1,2>(0,0);
        
        V << A-C,B-C;
        X_ed.block<2,2>(0,2*t) = V.cast<double>();
    }

    // updateDeformationMesh
    deformedState_st = new DeformationGradient2d(X_st * Bs);
    deformedState_ed = new DeformationGradient2d(X_ed * Bs);

    deformedState_st->svd();
    deformedState_ed->svd();
}



void TriangleMesh::addDeformationState(float t)
{
    cout << "====== add deformation state =====" << endl;

//    DeformationGradient2d newDeformationState = DeformationGradient2d();


    // interpolate in log space
    Eigen::MatrixXd U = (t * deformedState_st->_U.log() + (1 - t) * deformedState_ed->_U.log()).exp();
    Eigen::MatrixXd V = (t * deformedState_st->_V.log() + (1 - t) * deformedState_ed->_V.log()).exp();

    // interpolate singular
    Eigen::VectorXd s = (t * deformedState_st->_s.log() + (1 - t) * deformedState_ed->_s.log()).exp();

//    newDeformationState.updateRT(U, V, s);

    Eigen::MatrixXd S = s.asDiagonal();
    Eigen::MatrixXd Vt = V.transpose().conjugate();
    Eigen::MatrixXd R = U * Vt;
    Eigen::MatrixXd T = V * S * Vt;

    // compute the new mesh position
    int numV = initialVertices.rows();
    cout << "interpolate in log space: " << (1 - t) * deformedState_ed->_U.log() << endl;

    Eigen::MatrixXd deformedMesh = Eigen::MatrixXd(numV, 3);

//    deformedMesh = initialVertices * newDeformationState.R * newDeformationState.T;
    deformedMesh = initialVertices * R * T;

    // store new position in the vector
    predictedVertices.push_back(deformedMesh);

}


