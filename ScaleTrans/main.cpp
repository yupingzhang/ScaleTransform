/*
Input: rest mesh, and two deformation with the same topology

 1. create triangle mesh object, initialize rest position and deformed position

 2. compute deformation gradients, where dx, dy are the differences
 between the part position and the reference part location

 3. do SVD and get U, S, V

 4. compute log space interpolation for U, S, V (at t)

 5. compute R, T at t (rotation and stretch)

 6. compute new mesh position

Dependencies: Eigen and libigl

*/

#include <igl/file_exists.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>

#include <algorithm>
#include <iostream>
#include <unistd.h>

#include "TriangleMesh.h"


using namespace std;

int main(int argc, char *argv[])
{
    // if (argc < 2)
    // {
    //     cout << "Usage: ./xxx interp" << endl;
    // }
    // // interpolation in between
    // float p = stof(argv[1]);

    Eigen::MatrixXd V;
    Eigen::MatrixXd V1;
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F;

    // Load mesh
    cout << "load meshes..." << endl;
    if (!igl::file_exists("t0.obj") || !igl::file_exists("t_long.obj") || !igl::file_exists("t_tall.obj") ) {
        cout << "One or more files are not exist." << endl;
    }

    igl::read_triangle_mesh("t0.obj", V, F);
    igl::read_triangle_mesh("t_long.obj", V1, F);
    igl::read_triangle_mesh("t_tall.obj", V2, F);

    cout << "Base mesh: vert: " << V.rows() << " tri: " << F.rows() << endl;

    // create the mesh
    TriangleMesh* mesh = new TriangleMesh(V, V1, V2, F);
    mesh->initMesh();

    // add a new state to the mesh
    int numV = V.rows(), numT = F.rows();
    Eigen::MatrixXd* deformMesh = new Eigen::MatrixXd(numV, 3);
    Eigen::MatrixXd* deformGrad = new Eigen::MatrixXd(numT * 2, 2);

    for (int i = 0; i <= 10; ++i)
    {
        float p = float(i)/10.0;
        cout << "p: " << p << endl;
        
        mesh->addDeformationState(p, deformGrad);
        mesh->recoverMesh(0, deformGrad, deformMesh);
        
        // Create a libigl Viewer object
        igl::opengl::glfw::Viewer viewer;
        viewer.data().set_mesh(V1, F);
    
        // Set the vertices and faces for the viewer
        viewer.data().clear();
        viewer.data().set_mesh(*deformMesh, F);
        // Launch a viewer instance
        viewer.launch();
    }

    delete deformMesh;
    delete deformGrad;
   
    return 0;

}

