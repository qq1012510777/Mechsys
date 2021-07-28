#pragma once
#include "Intersection_Frac_boost.h"
#include <Eigen/Dense>
#include <deque>
#include <iostream>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Intersection_between_polygon
{
public:
    std::vector<Vector2d> Intersection;

public:
    Intersection_between_polygon(std::vector<Vector2d> F1, std::vector<Vector2d> F2);
};

inline Intersection_between_polygon::Intersection_between_polygon(std::vector<Vector2d> F1, std::vector<Vector2d> F2)
{
    std::vector<point_type> points_1(F1.size() + 1);

    for (size_t i = 0; i < F1.size(); ++i)
        points_1[i] = point_type(F1[i](0), F1[i](1));

    points_1[points_1.size() - 1] = point_type(F1[0](0), F1[0](1));

    polygon_type poly_1;
    boost::geometry::assign_points(poly_1, points_1);

    std::vector<point_type> points_2(F2.size() + 1);

    for (size_t i = 0; i < F2.size(); ++i)
        points_2[i] = point_type(F2[i](0), F2[i](1));

    points_2[points_2.size() - 1] = point_type(F2[0](0), F2[0](1));

    polygon_type poly_2;
    boost::geometry::assign_points(poly_2, points_2);

    std::deque<polygon_type> output;
    boost::geometry::correct(poly_1);
    boost::geometry::correct(poly_2);
    
    boost::geometry::intersection(poly_1, poly_2, output);

    if (output.size() > 1)
    {
        string ssy = "In class 'Intersection_between_polygon', intersection between two 2D polygons must be only one object!\nbut the size of vector of objects is: ";
        ssy = ssy + to_string(output.size());
        ssy = ssy + "\n";
        throw Error_throw_ignore(ssy);
    }

    if (output.size() == 0)
    {
        Intersection.clear();
        return;
    }

    for (auto it = boost::begin(boost::geometry::exterior_ring(output[0])); it != boost::end(boost::geometry::exterior_ring(output[0])); ++it)
    {
        double x = boost::geometry::get<0>(*it);
        double y = boost::geometry::get<1>(*it);

        Intersection.push_back(Vector2d{x, y});
    }

    Intersection.erase(Intersection.begin() + Intersection.size());
};
}; // namespace DFN