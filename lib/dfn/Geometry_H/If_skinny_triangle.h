#pragma once
#include "Eigen/Dense"
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{
class If_skinny_triangle
{
public:
    /* data */
    double angle = 0;
    bool If_skinny = false;

public:
    If_skinny_triangle(std::vector<RowVector3f> coord);
};

inline If_skinny_triangle::If_skinny_triangle(std::vector<RowVector3f> coord)
{
    RowVector3f V1 = coord[1] - coord[0],
                V2 = coord[2] - coord[0];

    double x1, x2, y1, y2, z1, z2;
    x1 = V1[0];
    y1 = V1[1];
    z1 = V1[2];

    x2 = V2[0];
    y2 = V2[1];
    z2 = V2[2];

    double dot = x1 * x2 + y1 * y2 + z1 * z2; //    #between [x1, y1, z1] and [x2, y2, z2]
    double lenSq1 = x1 * x1 + y1 * y1 + z1 * z1;
    double lenSq2 = x2 * x2 + y2 * y2 + z2 * z2;      //
    angle = acos(dot / sqrt(lenSq1 * lenSq2)); //

    angle = angle * 180 / M_PI;

    double E = 0.5;
    if (abs(angle - 0) < E ||
        abs(angle - 180) < E ||
        abs(angle - 360) < E)
        If_skinny = true;

    //cout << "angle: " << angle << ",  ";
    //cout << "abs(angle - 180): " << If_skinny << endl;
}

}; // namespace DFN