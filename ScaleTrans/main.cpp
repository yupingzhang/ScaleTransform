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
#include <igl/png/writePNG.h>

#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <cmath>

#include "TriangleMesh.h"
#define PI 3.1415926

using namespace std;


int idx;
float t;
Eigen::MatrixXd V;
Eigen::MatrixXi F;


void createDeformMesh0()
{
    Eigen::MatrixXd V;
    Eigen::MatrixXd V1;
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F;

    // Load mesh
    cout << "load meshes..." << endl;
    if (!igl::file_exists("t0.obj") || !igl::file_exists("t1.obj") || !igl::file_exists("t2.obj") ) {
        cout << "One or more files are not exist." << endl;
    }

    igl::read_triangle_mesh("t0.obj", V, F);
    igl::read_triangle_mesh("t1.obj", V1, F);
    igl::read_triangle_mesh("t2.obj", V2, F);

    cout << "Base mesh: vert: " << V.rows() << " tri: " << F.rows() << endl;

    // create the mesh
    TriangleMesh* mesh = new TriangleMesh(V, V1, V2, F);
    mesh->initMesh();

    // add a new state to the mesh
    int numV = V.rows(), numT = F.rows();
    Eigen::MatrixXd* deformMesh = new Eigen::MatrixXd(numV, 3);
    Eigen::MatrixXd* deformGrad = new Eigen::MatrixXd(numT * 2, 2);

    // for (int i = 0; i <= 10; ++i)
    // {
    int i = 8;
    float p = float(i)/10.0;
    cout << "p: " << p << endl;

    mesh->addDeformationState2d(p, deformGrad);
    mesh->recoverMesh2d(i, deformGrad, deformMesh);
    
    // Create a libigl Viewer object
    igl::opengl::glfw::Viewer viewer;
    // viewer.data().set_mesh(V1, F);

    // Set the vertices and faces for the viewer
    viewer.data().clear();
    viewer.data().set_mesh(*deformMesh, F);
    // Launch a viewer instance
    viewer.launch();
    // }

    delete deformMesh;
    delete deformGrad;

}


void createDeformMesh2d()
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    igl::read_triangle_mesh("t0.obj", V, F);

    TriangleMesh* mesh = new TriangleMesh(V, F);
    mesh->initMeshDeform2d();

    int numV = V.rows(), numT = F.rows();
    Eigen::MatrixXd* deformMesh = new Eigen::MatrixXd(numV, 3);
    Eigen::MatrixXd* deformGrad = new Eigen::MatrixXd(numT * 2, 2);

    for (int i = 0; i <= 10; ++i)
    {
        float p = float(i)/10.0;
        cout << "p: " << p << endl;

        mesh->addDeformationState2d(p, deformGrad);
        mesh->recoverMesh2d(i, deformGrad, deformMesh);
        
        // Create a libigl Viewer object
        igl::opengl::glfw::Viewer viewer;
        viewer.data().clear();
        viewer.data().set_mesh(*deformMesh, F);
        // Launch a viewer instance
        viewer.launch();
    }

    delete deformMesh;
    delete deformGrad;

}


bool pre_draw(igl::opengl::glfw::Viewer & viewer)
{
    viewer.data().set_mesh(V, F);
    viewer.core.align_camera_center(V, F);
    if (viewer.core.is_animating)
    {
        // Clear should be called before drawing the mesh
        viewer.data().clear();

        // Allocate temporary buffers for 1280x800 image
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

        int numV = V.rows(), numT = F.rows();
        Eigen::MatrixXd* deformMesh = new Eigen::MatrixXd(numV, 3);
        Eigen::MatrixXd* deformGrad = new Eigen::MatrixXd(numT * 3, 3);
           
        float p = float(idx++)/t;
        
        if (idx > t)
        {
            viewer.core.is_animating = false;
        }

        TriangleMesh* mesh = new TriangleMesh(V, F);
        mesh->initMeshDeform3d();
        mesh->addDeformationState3d(p, deformGrad);
        mesh->recoverMesh3d(idx, deformGrad, deformMesh);

        // camera_eye << 0.0, 0.0, 5.0
        float p_off = p + 0.2;
        viewer.core.camera_eye << 5 * cos(p_off * PI), 0.0, 5 * sin(p_off * PI);
        // light_position<< 0.0f, -0.30f, -5.0f;
        viewer.core.light_position << -5 * cos(p * PI), -0.30f, -5 * sin(p * PI);

        viewer.data().clear();
        viewer.data().set_mesh(*deformMesh, F);
        viewer.data().show_lines = false;
        // viewer.draw();

        // Draw the scene in the buffers
        viewer.core.draw_buffer(viewer.data(), true, R, G, B, A);

        // Save it to a PNG
        std::stringstream ss;
        ss << std::setw(3) << std::setfill('0') << idx;
        std::string s = ss.str();
        string filename = "out/out_" + s + ".png";
        igl::png::writePNG(R, G, B, A, filename);

        viewer.draw();

        delete deformMesh;
        delete deformGrad;
    }

    return false;
}


void createDeformMesh3d()
{
    cout << "createDeformMesh3d: load bunny.obj" << endl;
    igl::read_triangle_mesh("bunny.obj", V, F);

    cout << "Create a libigl Viewer object" << endl;
    igl::opengl::glfw::Viewer viewer;
    // viewer.callback_key_down = &key_down;
    viewer.callback_pre_draw = &pre_draw;
    viewer.core.is_animating = false;
    viewer.core.animation_max_fps = 1.;
    viewer.launch();
}


int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "Usage: ./xxx interp_total" << endl;
    }
    // interpolation in between
    t = stof(argv[1]);
    cout << "Total frame: " << t << endl;

    if (t > 100) 
    {
        return 0;
    }

    idx = 0;
    createDeformMesh3d();
}