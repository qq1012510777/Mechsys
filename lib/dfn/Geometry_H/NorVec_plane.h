#pragma once

#include "../Quaternion_H/Quaternion.h"
#include "Dense"
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{

//the function below finds the vector that
//(1) is vertical to fracture normal vector;
//and (2) lies on the horizontal plane (z = 0)

class NorVec_plane
{
public:
    Vector3d Normal_vector;
    NorVec_plane(const std::vector<Vector3d> _Verts);
};

inline NorVec_plane::NorVec_plane(const std::vector<Vector3d> _Verts)
{
    if (_Verts.size() < 3)
    {
        throw Error_throw_ignore("cannot calculate the normal vector of a plane with less than three point\n");
    }

    Vector3d Ns = (_Verts[0] - _Verts[1]).cross(_Verts[1] - _Verts[2]);

    if (Ns(2) < 0)
    {
        Ns = -Ns;
    };

    Normal_vector = Ns;
};
}; // namespace DFN