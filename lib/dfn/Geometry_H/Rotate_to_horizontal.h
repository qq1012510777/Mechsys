#pragma once
#include "../Geometry_H/Vector_2.h"
#include "../Quaternion_H/Quaternion.h"
#include "Eigen/Dense"
#include "Polygon_convex_3D.h"
#include "Rotation_verts.h"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Rotate_to_horizontal
{
public:
    Vector3d Ref_axis;
    double theta_PI;
    double theta_degree;
    double z;

public:
    // verts should be coplanar with Poly
    Rotate_to_horizontal(DFN::Polygon_convex_3D Poly);
    Rotate_to_horizontal(DFN::Polygon_convex_3D Poly, std::vector<Vector3d> &verts_rotate);
    void Rotate_other_pnts(const std::vector<Vector3d> verts, std::vector<Vector3d> &verts_rotate, bool &If_on_horizontal);
    void Rotate_back(const std::vector<Vector3d> verts, std::vector<Vector3d> &verts_rotate);
    void Rotate_back_without_z(const std::vector<Vector3d> verts, std::vector<Vector3d> &verts_rotate);
};

inline Rotate_to_horizontal::Rotate_to_horizontal(DFN::Polygon_convex_3D Poly)
{
    DFN::Vector_2 v(Poly.Normal_vector, Ref_axis);
    theta_degree = Poly.Beta;
    theta_PI = theta_degree * M_PI / 180;
}

inline Rotate_to_horizontal::Rotate_to_horizontal(DFN::Polygon_convex_3D Poly, std::vector<Vector3d> &verts_rotate)
{
    std::vector<Vector3d> verts = Poly.Corners;
    std::vector<Vector3d> normal_1(1), normal_2(1);
    normal_1[0] = Poly.Normal_vector;

    DFN::Vector_2 v(Poly.Normal_vector, Ref_axis);

    theta_degree = Poly.Beta;
    theta_PI = theta_degree * M_PI / 180;
    Quaternion_t Q_axis_1;

    if (theta_degree > 0.0001)
    {
        //cout << "1;\n";
        verts_rotate.resize(verts.size());
        DFN::Rotation_verts R1(verts, theta_PI, Q_axis_1, Ref_axis, verts_rotate);
        DFN::Rotation_verts R2(normal_1, theta_PI, Q_axis_1, Ref_axis, normal_2);
    }
    else
    {
        verts_rotate = verts;
        normal_2 = normal_1;
    }

    //-------a piece of debug code----
    for (size_t i = 0; i < verts_rotate.size(); i++)
    {
        if (abs(verts_rotate[i](2) - verts_rotate[(i + 1) % verts_rotate.size()](2)) > 0.0001)
        {
            DFN::Polygon_convex_3D ty1(verts_rotate);

            cout << "\npolygon 1\n";
            ty1.Print_Polygon();
            throw Error_throw_ignore("Error!!! The Z values of all vertexes should be the same! In class 'Rotate_to_horizontal'!\n");
        }
    }
    //--------------------------------
    this->z = verts_rotate[0](2);
    for (size_t i = 0; i < verts_rotate.size(); i++)
        verts_rotate[i](2) = 0;
}

inline void Rotate_to_horizontal::Rotate_other_pnts(const std::vector<Vector3d> verts, std::vector<Vector3d> &verts_rotate, bool &If_on_horizontal)
{
    Quaternion_t Q_axis_1;

    if (theta_degree > 0.0001)
    {
        //cout << "1;\n";
        verts_rotate.resize(verts.size());
        DFN::Rotation_verts R1(verts, theta_PI, Q_axis_1, Ref_axis, verts_rotate);
    }
    else
    {
        verts_rotate = verts;
    }

    for (size_t i = 0; i < verts_rotate.size(); i++)
        verts_rotate[i](2) -= z;

    If_on_horizontal = true;
    for (size_t i = 0; i < verts_rotate.size(); i++)
    {
        if (abs(verts_rotate[i](2)) > 0.05)
        {
            If_on_horizontal = false;
            break;
        }
    }
};

inline void Rotate_to_horizontal::Rotate_back(const std::vector<Vector3d> verts, std::vector<Vector3d> &verts_rotate)
{
    std::vector<Vector3d> verts_1 = verts;
    for (size_t i = 0; i < verts_1.size(); i++)
        verts_1[i](2) += z;

    Quaternion_t Q_axis_1;

    if (theta_degree > 0.0001)
    {
        //cout << "1;\n";
        verts_rotate.resize(verts_1.size());
        DFN::Rotation_verts R1(verts_1, -theta_PI, Q_axis_1, Ref_axis, verts_rotate);
    }
    else
    {
        verts_rotate = verts_1;
    }
};

inline void Rotate_to_horizontal::Rotate_back_without_z(const std::vector<Vector3d> verts, std::vector<Vector3d> &verts_rotate)
{
    Quaternion_t Q_axis_1;

    if (theta_degree > 0.0001)
    {
        //cout << "1;\n";
        verts_rotate.resize(verts.size());
        DFN::Rotation_verts R1(verts, -theta_PI, Q_axis_1, Ref_axis, verts_rotate);
    }
    else
    {
        verts_rotate = verts;
    }
};

}; // namespace DFN