#pragma once
#include "Dense"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Polygon_convex_3D
{
public:
    std::vector<Vector3d> Corners;
    Vector3d Normal_vector;
    double Beta; //degrees
    Vector4d Plane_parameter;

public:
    Polygon_convex_3D(const std::vector<Vector3d> Verts_1);
    Polygon_convex_3D();
    void Create(Polygon_convex_3D A);
    void Print_Polygon();
};

inline Polygon_convex_3D::Polygon_convex_3D()
{
    ;
};

inline void Polygon_convex_3D::Create(Polygon_convex_3D A)
{
    this->Corners = A.Corners;
    this->Normal_vector = A.Normal_vector;
    this->Beta = A.Beta;
    this->Plane_parameter = A.Plane_parameter;
};

inline Polygon_convex_3D::Polygon_convex_3D(const std::vector<Vector3d> Verts_1)
{
    if (Verts_1.size() < 3)
    {
        cout << "Polygon_convex_3D::Polygon_convex_3D, cannot create a polygon with less than three vertices!\n";
        exit(0);
    }
    Corners = Verts_1;
    Normal_vector = (Corners[0] - Corners[1]).cross(Corners[0] - Corners[2]);
    if (Normal_vector(2) < 0)
    {
        Normal_vector = -Normal_vector;
    }

    Normal_vector = Normal_vector / Normal_vector.norm();
    Beta = acos(Normal_vector(2)) * 180. / M_PI;
    if (Beta < -0.0001 || Beta > 90.0001)
    {
        cout << "Polygon_convex_3D::Polygon_convex_3D, Beta of a normal vector is incorrect!\n";
        cout << "Beta = " << Beta << endl;
        exit(0);
    }

    Plane_parameter << Normal_vector(0), Normal_vector(1), Normal_vector(2),
        -Normal_vector.dot(Corners[0]);
};

void Polygon_convex_3D::Print_Polygon()
{
    for (size_t i = 0; i < Corners.size(); ++i)
    {
        cout << Corners[i].transpose() << endl;
    }
};

}; // namespace DFN