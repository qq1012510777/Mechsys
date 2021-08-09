#pragma once
#include "../Error_throw/Error_throw.h"
#include <Eigen/Dense>
#include <deque>
#include <iostream>

using namespace std;
using namespace Eigen;

namespace DFN
{
class Triangle_orientation
{
public:
    bool If_clockwise = false;

public:
    Triangle_orientation(std::vector<Vector2d> polygon);
};

inline Triangle_orientation::Triangle_orientation(std::vector<Vector2d> polygon)
{

    if (polygon.size() != 3)
    {
        throw Error_throw_pause("in class 'Triangle_orientation', not a triangle!\n");
    }

    double p1x = polygon[0][0], p2x = polygon[1][0], p3x = polygon[2][0];
    double p1y = polygon[0][1], p2y = polygon[1][1], p3y = polygon[2][1];

    double val = (p2y - p1y) * (p3x - p2x) - (p2x - p1x) * (p3y - p2y);

    if (val > 0)
        If_clockwise = true;
    else if (val < 0)
        If_clockwise = false;
    else if (val == 0)
    {
        throw Error_throw_pause("in class 'Triangle_orientation', triangular element is a line segment!\n");  
    }
};

}; // namespace DFN