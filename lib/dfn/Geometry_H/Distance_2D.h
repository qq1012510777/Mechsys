#pragma once

#include "../Quaternion_H/Quaternion.h"
#include "Dense"
#include "Line_seg_2D.h"
#include "Point_2D.h"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Distance
{
public:
    Distance();
    double Distance_between_point_lineseg(const Point_2D p, const Line_seg_2D line, bool &if_p_on_line);
};

inline Distance::Distance()
{
    ;
};

inline double Distance::Distance_between_point_lineseg(const Point_2D p, const Line_seg_2D line, bool &if_p_on_line)
{
    std::vector<Vector3d> Verts{Vector3d{line.Point[0](0), line.Point[0](1), 0}, Vector3d{line.Point[1](0), line.Point[1](1), 0}};

    DFN::Point_2D p_c{p.Coordinate};
    bool fe = p_c.If_lies_on_the_bounds_of_polygon(Verts);

    if_p_on_line = false;

    if (fe == true)
    {
        if_p_on_line = true;
        return 1e10;
    }

    Vector2d B, C;
    B = line.Point[1] - line.Point[0];
    C = p.Coordinate - line.Point[0];

    if (abs(B(0)) < 0.0001 && abs(B(1)) < 0.0001)
    {
        return (B - C).norm();
    }
    else
    {
        double angle = atan2(B(1), B(0));
        Vector3d axis_z;
        axis_z << 0, 0, -1;

        Quaternion_t Q_axis;

        std::vector<Vector3d> ER(2), EU(2);
        ER[0] << B(0), B(1), 0;
        ER[1] << C(0), C(1), 0;

        DFN::Rotation_verts roe{ER, angle, Q_axis, axis_z, EU};

        if ((EU[1](0) >= 0 && EU[1](0) <= EU[0](0)) ||
            (EU[1](0) < 0 && EU[1](0) >= EU[0](0)))
        {
            return abs(EU[1](1));
        }
        else
        {
            double dist1 = abs(EU[1].norm());
            double dist2 = abs((EU[0] - EU[1]).norm());
            return dist1 < dist2 ? dist1 : dist2;
        }
    }
};

}; // namespace DFN