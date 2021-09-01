#pragma once
#include "Eigen/Dense"
#include <iostream>
using namespace std;
using namespace Eigen;

namespace DFN
{
class GeoOri_to_NormVec
{
public:
    GeoOri_to_NormVec(double dip_direction,
                      double dip_angle,
                      Vector3d &a);
};

inline GeoOri_to_NormVec::GeoOri_to_NormVec(double dip_direction,
                                            double dip_angle,
                                            Vector3d &a)
{
    ///spherical system firstly
    double alpha = 0, beta = 0;
    beta = dip_angle;
    if (dip_direction >= 90)
        alpha = 450 - dip_direction;
    else if (dip_direction <= 90)
        alpha = 90 - dip_direction;

    //------------------
    a(0) = sin(beta * M_PI / 180) * cos(alpha / 180.0 * M_PI);
    a(1) = sin(beta / 180.0 * M_PI) * sin(alpha / 180.0 * M_PI);
    a(2) = cos(beta / 180.0 * M_PI);

    if (a(2) < 0)
    {
        a = -a;
    }
};

} // namespace DFN