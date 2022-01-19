#pragma once

#include "Eigen/Dense"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

//the function below finds the vector that
//(1) is vertical to fracture normal vector;
//and (2) lies on the horizontal plane (z = 0)

class Interval_1D
{
public:
    double x1_min;
    double x1_max;
    double x2_min;
    double x2_max;
    bool if_intersect;
    Vector2d Intersection;

public:
    Interval_1D(const Vector2d A, const Vector2d B);
    void Interval_1D_Known(Vector2d &intersection);
};

inline Interval_1D::Interval_1D(const Vector2d A, const Vector2d B)
{
    if (A(0) > A(1))
    {
        x1_min = A(1);
        x1_max = A(0);
    }
    else
    {
        x1_min = A(0);
        x1_max = A(1);
    }

    if (B(0) > B(1))
    {
        x2_min = B(1);
        x2_max = B(0);
    }
    else
    {
        x2_min = B(0);
        x2_max = B(1);
    }

    if (x1_min - x2_max > 0.0001 || //tolerance
        x2_min - x1_max > 0.0001)
    {
        if_intersect = false;
        return;
    }
    else
    {
        Vector2d Intersection_c;
        if (x2_max > x1_max)
        {
            Interval_1D_Known(Intersection_c);
        }
        else
        {
            Interval_1D_Known(Intersection_c);
        }
        Intersection = Intersection_c;
        if_intersect = true;
        return;
    }
};

inline void Interval_1D::Interval_1D_Known(Vector2d &intersection)
{
    //we know that B(1) >= A(1)
    Vector2d intersection_ts;
    intersection_ts(1) = x1_max < x2_max ? x1_max : x2_max;

    if (x2_min > x1_min)
    {

        intersection_ts(0) = x2_min;
    }
    else
    {
        intersection_ts(0) = x1_min;
    }
    
    double AF = intersection_ts(0) - intersection_ts(1);
    
    if (abs(AF) < 0.0001)
    {
        intersection(0) = intersection_ts(0);
        intersection(1) = intersection_ts(0);
    }
    else
        intersection = intersection_ts;
};

}; // namespace DFN