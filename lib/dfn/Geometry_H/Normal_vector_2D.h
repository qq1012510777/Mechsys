#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
using namespace Eigen;

namespace DFN
{
class Normal_vector_2D
{
public:
    Vector2d Normal;
    //double cos_theta_y;
    //double cos_theta_x;

public:
    Normal_vector_2D();
    void Normal_vector_of_a_line(std::vector<Vector2d> Line);
};

Normal_vector_2D::Normal_vector_2D(){
    //
};

void Normal_vector_2D::Normal_vector_of_a_line(std::vector<Vector2d> Line)
{
    Normal << -(Line[1][1] - Line[0][1]), Line[1][0] - Line[0][0];
    double norm = Normal.norm();
    Normal[0] = Normal[0] / norm;
    Normal[1] = Normal[1] / norm;

    //cos_theta_x = abs(Normal[0]) / Normal.norm();
    //cos_theta_y = abs(Normal[1]) / Normal.norm();
};
}; // namespace DFN