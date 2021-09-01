#pragma once
#include "Eigen/Dense"
#include <iostream>
using namespace std;
using namespace Eigen;

namespace DFN
{
class NormVec_to_GeoOri
{
public:
    NormVec_to_GeoOri(Vector3d a,
                      double &dip_direction,
                      double &dip_angle);
};

inline NormVec_to_GeoOri::NormVec_to_GeoOri(Vector3d a,
                                            double &dip_direction,
                                            double &dip_angle)
{
    double l = a[0], m = a[1], n = a[2];

    if (l > 1.0 || m > 1.0 || n > 1.0)
    {
        double max = l > m ? l : m;
        max = max > n ? max : n;

        l = l / max;
        m = m / max;
        n = n / max;
    }

    double beta_tmp = acos(n) * 180.0 / M_PI;
    double alpha_tmp = atan2(m, l) * 180.0 / M_PI;

    if (alpha_tmp < 0)
        alpha_tmp = 360 + alpha_tmp;

    dip_angle = beta_tmp;

    if (alpha_tmp <= 90)
        dip_direction = 90 - alpha_tmp;
    else if (alpha_tmp > 90)
        dip_direction = 450 - alpha_tmp;
};

} // namespace DFN