#pragma once

#include "../Quaternion_H/Quaternion.h"
#include "Dense"
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

class Vector_2
{

public:
    Vector_2(const Vector3d Normal_vector, Vector3d &temp3);
};

inline Vector_2::Vector_2(const Vector3d Normal_vector, Vector3d &temp3)
{
    Vector3d project_xy;
    project_xy << Normal_vector(0), Normal_vector(1), 0.;

    Vector3d axis;
    axis << 0, 0, -1.;

    double R_angle_temp1 = 90.0 * M_PI / 180.0;

    Quaternion_t Q_axis;
    NormalizeRotation(R_angle_temp1, axis, Q_axis);
    Rotation(project_xy, Q_axis, temp3);
};

}; // namespace DFN