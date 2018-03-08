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
    int numTri = triangles.rows();    // get triangle number
    Bs.resize(2*numTri, 2); 
    F_st.resize(2*numTri, 2);
    F_ed.resize(2*numTri, 2);

    Eigen::Vector2d A,B,C;
    
    for(int t=0; t<numTri; t++)
    {
        Eigen::Matrix2d dx1, dx2, X0, beta;

        A = initialVertices.row(triangles(t,0)).block<1,2>(0,0);
        B = initialVertices.row(triangles(t,1)).block<1,2>(0,0);
        C = initialVertices.row(triangles(t,2)).block<1,2>(0,0);

        X0 << A-C,B-C;
        beta = X0.inverse();
        Bs.block<2,2>(2*t, 0) = X0;

        ////// deformation start
        A = deformedVertices_st.row(triangles(t,0)).block<1,2>(0,0);
        B = deformedVertices_st.row(triangles(t,1)).block<1,2>(0,0);
        C = deformedVertices_st.row(triangles(t,2)).block<1,2>(0,0);

        dx1 << A-C,B-C; 
        F_st.block<2,2>(2*t, 0) = dx1 * beta;

        ////// deformation end
        A = deformedVertices_ed.row(triangles(t,0)).block<1,2>(0,0);
        B = deformedVertices_ed.row(triangles(t,1)).block<1,2>(0,0);
        C = deformedVertices_ed.row(triangles(t,2)).block<1,2>(0,0);

        dx2 << A-C,B-C;
        F_ed.block<2,2>(2*t, 0) = dx2 * beta;
    }

    // updateDeformationMesh
    deformedState_st = new DeformationGradient2d(numTri, F_st);
    deformedState_ed = new DeformationGradient2d(numTri, F_ed);

    deformedState_st->svd();
    deformedState_ed->svd();
}



void TriangleMesh::addDeformationState(float p, Eigen::MatrixXd* deformation)
{
    cout << "====== add deformation state =====" << endl;
    Eigen::MatrixXd U = Eigen::Matrix2d::Zero();
    Eigen::MatrixXd S = Eigen::Matrix2d::Zero();
    Eigen::MatrixXd Vt = Eigen::Matrix2d::Zero();

    int numTri = triangles.rows();
    for (int t = 0; t < numTri; ++t)
    {
        const double u = p * deformedState_st->theta_u[t] + (1 - p) * deformedState_ed->theta_u[t];
        const double v = p * deformedState_st->theta_v[t] + (1 - p) * deformedState_ed->theta_v[t];

        Eigen::Rotation2Dd r_u(u); 
        Eigen::Rotation2Dd r_v(v); 

        Eigen::Matrix2d U = Eigen::Matrix2d(r_u.toRotationMatrix());
        Eigen::Matrix2d V = Eigen::Matrix2d(r_v.toRotationMatrix());

        Eigen::Matrix2d Vt = V.transpose().conjugate();

        // the stretch matrix interpolated in log space
        S(0, 0) = exp(p * log(deformedState_st->_s(t, 0)) + (1 - p) * log(deformedState_ed->_s(t, 0)));
        S(1, 1) = exp(p * log(deformedState_st->_s(t, 1)) + (1 - p) * log(deformedState_ed->_s(t, 1)));
        // 
        // S(0, 0) = p * deformedState_st->_s(t, 0) + (1-p) * deformedState_ed->_s(t, 0);
        // S(1, 1) = p * deformedState_st->_s(t, 1) + (1-p) * deformedState_ed->_s(t, 1);

        deformation->block<2,2>(2*t, 0) = U * S * Vt;
    }
}

void TriangleMesh::recoverMesh(int i, Eigen::MatrixXd* deformGrad, Eigen::MatrixXd* deformMesh)
{
    int numTri = triangles.rows();

    for(int t=0; t<numTri; t++) 
    {
        Eigen::Vector2d xy1 = deformGrad->block<2, 2>(2*t, 0) * initialVertices.row(triangles(t, 0)).transpose();
        Eigen::Vector2d xy2 = deformGrad->block<2, 2>(2*t, 0) * initialVertices.row(triangles(t, 1)).transpose();
        Eigen::Vector2d xy3 = deformGrad->block<2, 2>(2*t, 0) * initialVertices.row(triangles(t, 2)).transpose();
        
        deformMesh->row(triangles(t, 0)) = Eigen::Vector3d(xy1(0), xy1(1), 0.0);
        deformMesh->row(triangles(t, 1)) = Eigen::Vector3d(xy2(0), xy2(1), 0.0);
        deformMesh->row(triangles(t, 2)) = Eigen::Vector3d(xy3(0), xy3(1), 0.0);
    }
    
    cout << "==============  write to out.obj  ==============" << endl << deformMesh->rows() << " x " << deformMesh->cols() << endl;
    igl::writeOBJ("out/out_" + to_string(i) + ".obj", *deformMesh, triangles);

}
