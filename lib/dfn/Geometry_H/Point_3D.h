#pragma once
#include "Point_2D.h"
#include "Polygon_convex_3D.h"
#include "Rotate_to_horizontal.h"

typedef boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian> point_type_3D;
typedef boost::geometry::model::linestring<point_type_3D> Linestring_3D;
typedef boost::geometry::model::polygon<point_type_3D> polygon_type_3D;

namespace DFN
{

class Point_3D
{
public:
    Vector3d Coordinate;

public:
    Point_3D();
    Point_3D(Vector3d A);
    void Re_constructor(Vector3d A);

    bool If_lies_on_a_line_seg(const std::vector<Vector3d> Line_seg);
    bool If_lies_on_the_bounds_of_polygon(const vector<Vector3d> Verts_1);
    bool If_lies_within_a_polygon_3D(const vector<Vector3d> Verts_1);
    void Remove_the_same_pnt(std::vector<DFN::Point_3D> &PntSet);
    bool If_the_same_pnt(DFN::Point_3D S);
};

inline Point_3D::Point_3D(){};

inline Point_3D::Point_3D(Vector3d A)
{
    Coordinate = A;
};

inline void Point_3D::Re_constructor(Vector3d A)
{
    Coordinate = A;
};

inline bool Point_3D::If_lies_on_a_line_seg(const std::vector<Vector3d> Line_seg)
{
    if (Line_seg.size() != 2)
    {
        throw Error_throw_ignore("Input should be a line! In class 'Point_3D', function 'If_lies_on_a_line_seg'\n");
    }

    if ((Line_seg[0] - Line_seg[1]).norm() < 0.01)
    {
        if ((Coordinate - Line_seg[0]).norm() < 0.01)
        {
            return true;
        }
        else
            return false;
    }
    else
    {
        Linestring_3D EDGE_POLYGON{{Line_seg[0](0), Line_seg[0](1), Line_seg[0](2)},
                                   {Line_seg[1](0), Line_seg[1](1), Line_seg[1](2)}};

        point_type_3D p(Coordinate(0), Coordinate(1), Coordinate(2));

        double distanceA = boost::geometry::distance(p, EDGE_POLYGON);

        if (distanceA < 0.01)
            return true;

        return false;
    }
};

inline bool Point_3D::If_lies_on_the_bounds_of_polygon(const vector<Vector3d> Verts_1)
{
    for (size_t i = 0; i < Verts_1.size(); ++i)
    {
        std::vector<Vector3d> Lineseg = {Verts_1[i], Verts_1[(i + 1) % Verts_1.size()]};

        if (If_lies_on_a_line_seg(Lineseg) == true)
            return true;
    }

    return false;
};

inline bool Point_3D::If_lies_within_a_polygon_3D(const vector<Vector3d> Verts_1)
{
    DFN::Polygon_convex_3D poly{Verts_1};

    vector<Vector3d> Verts_2;
    DFN::Rotate_to_horizontal R1{poly, Verts_2};

    vector<Vector3d> Verts_3;
    bool On_H;
    R1.Rotate_other_pnts(vector<Vector3d>{this->Coordinate}, Verts_3, On_H);

    if (On_H == false)
        return false;
    else
    {
        Point_2D AS{Verts_3[0](0), Verts_3[0](1)};

        std::vector<Vector2d> Verts_4(Verts_2.size());
        for (size_t i = 0; i < Verts_4.size(); ++i)
            Verts_4[i] << Verts_2[i](0), Verts_2[i](1);

        if (AS.If_lies_within_a_polygon_boost(Verts_4) == true)
        {
            return true;
        }
        else
            return false;
    }
};

inline bool Point_3D::If_the_same_pnt(DFN::Point_3D S)
{
    if ((this->Coordinate - S.Coordinate).norm() < 0.05)
        return true;
    else
        return false;
};

inline void Point_3D::Remove_the_same_pnt(std::vector<DFN::Point_3D> &PntSet)
{
    std::vector<size_t> Index(PntSet.size());
    for (size_t i = 0; i < PntSet.size(); ++i)
        Index[i] = 0;

    for (size_t i = 0; i < PntSet.size() - 1; ++i)
    {
        if (Index[i] == 0)
        {
            for (size_t j = i + 1; j < PntSet.size(); ++j)
            {
                if (Index[j] == 0)
                {
                    if (PntSet[i].If_the_same_pnt(PntSet[j]) == true)
                    {
                        Index[j] = 1;
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < Index.size();)
    {
        if (Index[i] == 1)
        {
            Index.erase(Index.begin() + i);
            PntSet.erase(PntSet.begin() + i);
        }
        else
            ++i;
    }
};

} // namespace DFN