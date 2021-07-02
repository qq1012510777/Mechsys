#pragma once
#include "../Quaternion_H/Quaternion.h"
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

class Rotation_verts
{

public:
    Rotation_verts(vector<Vector3d> f1,
                   double R_angle_temp1,
                   Quaternion_t &Q_axis_1,
                   Vector3d temp1,
                   std::vector<Vector3d> &Verts_1);
};

inline Rotation_verts::Rotation_verts(vector<Vector3d> f1,
                                      double R_angle_temp1,
                                      Quaternion_t &Q_axis_1,
                                      Vector3d temp1,
                                      std::vector<Vector3d> &Verts_1)
{
    NormalizeRotation(R_angle_temp1, temp1, Q_axis_1);
    for (size_t i = 0; i < f1.size(); ++i)
    {
        Rotation(f1[i], Q_axis_1, Verts_1[i]);
    }
};

}; // namespace DFN