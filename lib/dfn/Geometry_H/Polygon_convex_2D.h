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

class Polygon_convex_2D
{
public:
    std::vector<Vector2d> Corners;
    double x_min;
    double x_max;
    double y_min;
    double y_max;

public:
    Polygon_convex_2D(const std::vector<Vector3d> Verts_1);
    void Extream();
};

inline Polygon_convex_2D::Polygon_convex_2D(const std::vector<Vector3d> Verts_1)
{
    Corners.resize(Verts_1.size());
    for (size_t i = 0; i < Verts_1.size(); ++i)
        Corners[i] << Verts_1[i](0), Verts_1[i](1);
};

inline void Polygon_convex_2D::Extream()
{
    x_min = Corners[0](0);
    x_max = Corners[0](0);

    y_min = Corners[0](1);
    y_max = Corners[0](1);
    for (size_t i = 1; i < Corners.size(); ++i)
    {
        if (Corners[i](0) < x_min)
            x_min = Corners[i](0);

        if (Corners[i](0) > x_max)
            x_max = Corners[i](0);

        if (Corners[i](1) < y_min)
            y_min = Corners[i](1);

        if (Corners[i](1) > y_max)
            y_max = Corners[i](1);
    }
};

}; // namespace DFN