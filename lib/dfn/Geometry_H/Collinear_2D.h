#pragma once

#include "Dense"
#include "Line_seg_2D.h"
#include "Rotation_verts.h"
#include "Vector_2.h"
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



class Collinear_2D
{
public:
    Collinear_2D(const Vector2d A,
                 const Vector2d B,
                 const Vector2d C,
                 bool &collinear);
};

inline Collinear_2D::Collinear_2D(const Vector2d A,
                                  const Vector2d B,
                                  const Vector2d C,
                                  bool &collinear)
{
    collinear = false;
    if (abs(A(0) - B(0)) < 0.0001 &&
        abs(A(0) - C(0)) < 0.0001)
    {
        collinear = true;
        return;
    }

    if (abs(A(1) - B(1)) < 0.0001 &&
        abs(A(1) - C(1)) < 0.0001)
    {
        collinear = true;
        return;
    }

    Vector2d YU;
    YU = B - A;

    if (abs(YU(0)) < 0.0001 && abs(YU(1)) < 0.0001)
    {
        //actually a point
        Vector2d IU;
        IU = A - C;
        if (abs(IU(0)) < 0.0001 && abs(IU(1)) < 0.0001)
        {
            collinear = true;
        }
        else
            collinear = false;
    }
    else
    {
        // rotataion
        double angle_radian = atan2(YU(1), YU(0));

        Vector3d axis{0, 0, -1};
        Vector3d C1{C(0) - A(0), C(1) - A(1), 0};
        std::vector<Vector3d> RF{C1};
        std::vector<Vector3d> RF2(1);

        Quaternion_t Q_axis;
        //cout << angle_radian * 180. / M_PI << endl;
        if (abs(YU(0)) < 0.0001)
        {
            angle_radian = 90.0 * M_PI / 180.0;
        }

        if (abs(YU(1)) > 0.0001)
        {
            //cout << "12" << endl;
            DFN::Rotation_verts R1(RF, angle_radian, Q_axis, axis, RF2);
        }
        else
            RF2 = RF;
        //cout << "RF2[0](1): " << RF2[0](1) << endl;
        if (abs(RF2[0](1)) < 0.0001)
        {
            collinear = true;
        }
        else
        {
            collinear = false;
        }
    }
    return;
};

}; // namespace DFN