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

class Line_seg_3D
{
public:
    vector<Vector3d> Point;
    Vector3d directional_vector;

public:
    Line_seg_3D(const Vector3d P1, const Vector3d P2);
    Vector3d z_Coordin_of_a_Point_in_the_Line_seg_3D(const double z = 0);
    bool Line_crosses_a_horizontal_plane(Vector3d &A, const double z = 0);
};

inline Line_seg_3D::Line_seg_3D(const Vector3d P1, const Vector3d P2)
{
    Point.resize(2);
    Point[0] = P1;
    Point[1] = P2;
    directional_vector = P2 - P1;
};

inline Vector3d Line_seg_3D::z_Coordin_of_a_Point_in_the_Line_seg_3D(const double z)
{
    //before this cuntion, you have to make sure that
    //the point must be belonging to the line.
    double x = 0, y = 0;

    if (abs(directional_vector(0)) > 0.0001 && abs(directional_vector(1)) > 0.0001)
    {
        x = (z - Point[0](2)) / directional_vector(2) * directional_vector(0) + Point[0](0);
        y = (z - Point[0](2)) / directional_vector(2) * directional_vector(1) + Point[0](1);
    }
    else if (abs(directional_vector(0)) < 0.0001 && abs(directional_vector(1)) > 0.0001)
    {
        //the line is parallel to y axis
        x = Point[0](0);
        y = (z - Point[0](2)) / directional_vector(2) * directional_vector(1) + Point[0](1);
    }
    else if (abs(directional_vector(0)) > 0.0001 && abs(directional_vector(1)) < 0.0001)
    {
        //the line is parallel to x axis
        x = (z - Point[0](2)) / directional_vector(2) * directional_vector(0) + Point[0](0);
        y = Point[0](1);
    }
    else if (abs(directional_vector(0)) < 0.0001 && abs(directional_vector(1)) < 0.0001)
    {
        x = Point[0](0); 
        y = Point[0](1);  
    }
    else
    {
        cout << "Line_seg_3D, Undefined behavior!\n";
        cout << directional_vector.transpose() << endl;
        cout << Point[0].transpose() << endl;
        cout << Point[1].transpose() << endl;
        exit(0);
    }

    Vector3d SO;
    SO << x, y, z;
    return SO;
};

inline bool Line_seg_3D::Line_crosses_a_horizontal_plane(Vector3d &A, const double z)
{
    double z_min = Point[0](2) <= Point[1](2)? Point[0](2):Point[1](2);
    double z_max = Point[0](2) >= Point[1](2)? Point[0](2):Point[1](2);
    z_min -= z;
    z_max -= z;

    if (z_min < 0 && z_max > 0)
    {
        A = z_Coordin_of_a_Point_in_the_Line_seg_3D(z);
        return true;
    }
    else
    {
        return false;
    }
};

}; // namespace DFN