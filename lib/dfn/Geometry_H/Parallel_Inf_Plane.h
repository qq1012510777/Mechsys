#pragma once
#include "../DFN_H/Fracture_WL.h"
#include "Dense"
#include "Polygon_convex_3D.h"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Parallel_Inf_Plane
{
public:
    string Raltion;

public:
    Parallel_Inf_Plane(const Polygon_convex_3D F1, const Polygon_convex_3D F2);
};

inline Parallel_Inf_Plane::Parallel_Inf_Plane(const Polygon_convex_3D F1, const Polygon_convex_3D F2)
{
    double a1 = F1.Plane_parameter(0);
    double b1 = F1.Plane_parameter(1);
    double c1 = F1.Plane_parameter(2);
    double d1 = F1.Plane_parameter(3);
    double a2 = F2.Plane_parameter(0);
    double b2 = F2.Plane_parameter(1);
    double c2 = F2.Plane_parameter(2);
    double d2 = F2.Plane_parameter(3);

    double index1 = 0; //when index1 = 1, means parallel; if index1 = 0, means two infinite planes intersect
    //double index2 = 0; //0 means parallel but not overlap, 1 means overlap

    Vector3d A, B;
    A << a1, b1, c1;
    B << a2, b2, c2;

    if (abs((A.cross(B)).norm()) < 0.0001)
    {
        index1 = 1;
        Raltion = "Parallel";
    }
    else
    {
        index1 = 0;
        Raltion = "Intersecting";
    }

    if (index1 == 1)
    {
        double x = 1.5, y = 3.6, z1 = 0, z2 = 0;
        z1 = -(a1 * x + b1 * y + d1) / c1;
        z2 = -(a2 * x + b2 * y + d2) / c2;
        if (abs(z1 - z2) < 0.0001)
        {
            //index2 = 1;
            Raltion = "Overlapping";
        }
    }
};
}; // namespace DFN