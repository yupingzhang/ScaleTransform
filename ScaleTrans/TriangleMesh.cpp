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
}

/**
 * { init the base matrix and two deformation gradient }
 */
void TriangleMesh::initMesh()
{
    cout << "====== initialize the mesh data =====" << endl;
    // Compute the base matrix for deformation gradients
    int numTris = triangles.rows();    // get triangle number
    Bs.resize(2*numTris, 2);
//    X0.resize(2, 2*numTris);
    X_st.resize(2*numTris, 2);
    X_ed.resize(2*numTris, 2);

    for(int t=0;t<numTris;t++)
    {
        Eigen::Vector2d A,B,C;
        Eigen::Matrix2d V, beta;

        A = initialVertices.row(triangles(t,0)).block<1,2>(0,0);
        B = initialVertices.row(triangles(t,1)).block<1,2>(0,0);
        C = initialVertices.row(triangles(t,2)).block<1,2>(0,0);

        beta << A-C,B-C;
        Bs.block<2,2>(2*t, 0) = beta.inverse();
//        X0.block<2,2>(0,2*t) = beta;

        ////// deformation start
        A = deformedVertices_st.row(triangles(t,0)).block<1,2>(0,0);
        B = deformedVertices_st.row(triangles(t,1)).block<1,2>(0,0);
        C = deformedVertices_st.row(triangles(t,2)).block<1,2>(0,0);

        V << A-C,B-C;
        X_st.block<2,2>(2*t, 0) = V * beta;

        ////// deformation end
        A = deformedVertices_ed.row(triangles(t,0)).block<1,2>(0,0);
        B = deformedVertices_ed.row(triangles(t,1)).block<1,2>(0,0);
        C = deformedVertices_ed.row(triangles(t,2)).block<1,2>(0,0);

        V << A-C,B-C;
        X_ed.block<2,2>(2*t, 0) = V * beta;
    }

    // updateDeformationMesh
    deformedState_st = new DeformationGradient2d(X_st);
    deformedState_ed = new DeformationGradient2d(X_ed);

    deformedState_st->svd();
    deformedState_ed->svd();
}



void TriangleMesh::addDeformationState(float t, Eigen::MatrixXd* deformation)
{
    cout << "====== add deformation state =====" << endl;

    // interpolate in log space
    Eigen::MatrixXd U = (t * deformedState_st->_U.log() + (1 - t) * deformedState_ed->_U.log()).exp();
    Eigen::MatrixXd V = (t * deformedState_st->_V.log() + (1 - t) * deformedState_ed->_V.log()).exp();

    // interpolate singular
    Eigen::VectorXd s = (t * deformedState_st->_s.log() + (1 - t) * deformedState_ed->_s.log()).exp();

    Eigen::MatrixXd S = s.asDiagonal();
    cout << "S: " << endl << S << endl;
    cout << "U rows: " << endl << U.rows() << " x " << U.cols() << endl;

    Eigen::MatrixXd Vt = V.transpose().conjugate();
    cout << "Vt: " << endl << Vt << endl;

    Eigen::MatrixXd R = U * Vt;              // ???
    Eigen::MatrixXd T = V * S * Vt;

    cout << "R: " << endl << R.rows() << " x " << R.cols() << endl;
    cout << "T: " <<endl << T.rows() << " x " << T.cols() << endl;

    *deformation = R * T;
    cout << "deformation: " << deformation->rows() << " x " << deformation->cols() << endl;

}

void TriangleMesh::recoverMesh(Eigen::MatrixXd* deformGrad, Eigen::MatrixXd* deformMesh)
{
//    Eigen::MatrixXd X = Eigen::MatrixXd(X0);

    int numTris = triangles.rows();
    for(int t=0; t<numTris; t++) {
        Eigen::Vector2d A, B, C;

        Eigen::Vector2d x0 = initialVertices.row(triangles(t, 0)).block<1, 2>(0, 0);
        Eigen::MatrixXd dx = deformGrad->row(t*2).block<2, 2>(0, 0);
        cout << dx << endl << endl;

//        A = x0;
//        B = A - Eigen::Vector2d(dx.row(1));
//        C = A - Eigen::Vector2d(dx.row(0));
        A = B = C = x0;

        //Todo: check if value is the same if overwrite
        deformMesh->row(triangles(t, 0)) = A;
        deformMesh->row(triangles(t, 1)) = B;
        deformMesh->row(triangles(t, 2)) = C;
    }

    cout << "==============  write to out.obj  ==============" << endl << deformMesh->rows() << " x " << deformMesh->cols() << endl;
    igl::writeOBJ("out.obj", *deformMesh, triangles);

}
