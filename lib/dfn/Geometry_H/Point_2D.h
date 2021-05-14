#pragma once
#include "Collinear_2D.h"
#include "Dense"
#include "Rotation_verts.h"
#include <cmath>
#include <iostream>
#include <list>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

#include "boost/assign/std/vector.hpp"
#include "boost/geometry.hpp"
#include "boost/geometry/algorithms/area.hpp"
#include "boost/geometry/algorithms/assign.hpp"
#include "boost/geometry/algorithms/intersection.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/geometries/polygon.hpp"
#include "boost/geometry/io/dsv/write.hpp"

using namespace boost::assign;
typedef boost::geometry::model::d2::point_xy<double> point_type;
typedef boost::geometry::model::polygon<point_type> polygon_type;
typedef boost::geometry::model::segment<point_type> segment_type;

namespace DFN
{

class Point_2D
{
public:
    Vector2d Coordinate;

public:
    Point_2D(const Vector2d A);
    Point_2D();
    void Re_constructor(const Vector2d A);
    bool If_lies_within_a_polygon_2D(const vector<Vector3d> Verts_1);
    bool If_lies_on_the_bounds_of_polygon(const vector<Vector3d> Verts_1);
    bool If_overlap_with_another_pnt(Point_2D AN);
    bool If_lies_on_a_line_seg(const std::vector<Vector2d> Line_seg);
    Vector2d Perpend_foot_on_a_line_seg(const std::vector<Vector2d> Line_seg);
};

inline Point_2D::Point_2D(const Vector2d A)
{
    Coordinate = A;
};

inline Point_2D::Point_2D()
{
    ;
};

inline void Point_2D::Re_constructor(const Vector2d A)
{
    Coordinate = A;
};

inline bool Point_2D::If_lies_within_a_polygon_2D(const vector<Vector3d> Verts_1)
{
    std::vector<point_type> points(Verts_1.size() + 1);
    for (size_t i = 0; i < Verts_1.size(); ++i)
        points[i] = point_type(Verts_1[i](0), Verts_1[i](1));

    points[points.size() - 1] = point_type(Verts_1[0](0), Verts_1[0](1));

    polygon_type poly;
    boost::geometry::assign_points(poly, points);

    point_type p(Coordinate(0), Coordinate(1));

    bool ty = boost::geometry::within(p, poly);

    return ty;
};

inline bool Point_2D::If_lies_on_the_bounds_of_polygon(const vector<Vector3d> Verts_1)
{
    bool EFF = false;
    Vector2d p;
    p << Coordinate(0), Coordinate(1);

    for (size_t i = 0; i < Verts_1.size(); ++i)
    {
        Vector2d A, B;
        A << Verts_1[i](0), Verts_1[i](1);
        B << Verts_1[(i + 1) % Verts_1.size()](0), Verts_1[(i + 1) % Verts_1.size()](1);

        if (abs((A - B)(0)) < 0.0001 && abs((A - B)(1)) < 0.0001)
        {
            Vector2d TF;
            TF = p - A;

            if (abs(TF(0)) < 0.0001 && abs(TF(1)) < 0.0001)
            {
                EFF = true;
                break;
            }
        }

        bool collinear1, collinear2;
        DFN::Collinear_2D stq1{p, A, B, collinear1};
        DFN::Collinear_2D stq2{A, B, p, collinear2};

        if (collinear1 == true && collinear2 == true)
        {
            Vector2d A2, B2;
            A2 = B - A;
            B2 = p - A;

            double angle_tmp = atan2(A2(1), A2(0));

            double angle_check = atan2(B2(1), B2(0));

            Quaternion_t Q_axis_1;
            Vector3d temp1;
            std::vector<Vector3d> Verts3(1), Verts4(1), Verts5(1), Verts6(1);
            temp1 << 0, 0, -1;

            Verts3[0] << A2(0), A2(1), 0;
            Verts5[0] << B2(0), B2(1), 0;

            if (round(abs(A2(1)), 4) > 0.0001)
            {
                DFN::Rotation_verts Re(Verts3, angle_tmp, Q_axis_1, temp1, Verts4);
                DFN::Rotation_verts Re2(Verts5, angle_tmp, Q_axis_1, temp1, Verts6);
            }
            else
            {
                Verts4 = Verts3;
                Verts6 = Verts5;
            }

            if (abs(Verts4[0](1)) > 0.0001 ||
                abs(Verts6[0](1)) > 0.0001)
            {
                if (abs((angle_tmp * 180. / M_PI) - (angle_check * 180. / M_PI)) > 1.0)
                {
                    double angle1, angle2;
                    if ((angle_tmp * 180. / M_PI) > 0)
                    {
                        angle1 = (angle_tmp * 180. / M_PI);
                    }
                    else
                    {
                        angle1 = (angle_tmp * 180. / M_PI) + 180.0;
                    }

                    if ((angle_check * 180. / M_PI) > 0)
                    {
                        angle2 = (angle_check * 180. / M_PI);
                    }
                    else
                    {
                        angle2 = (angle_check * 180. / M_PI) + 180.0;
                    }

                    if (abs(angle1 - angle2) > 1.0)
                    {
                        cout << "Point_2D::If_lies_on_the_bounds_of_polygon, Rotating to x_axis failed!\n";
                        cout << "Rotation angle: " << (angle_tmp * 180. / M_PI) << endl;
                        cout << "Rotation angle 2: " << (angle_check * 180. / M_PI) << endl;
                        cout << "\nThe three points, A, B and C are (before move):\n";
                        cout << A.transpose() << endl;
                        cout << B.transpose() << endl;
                        cout << p.transpose() << endl;

                        cout << "\nThe two points, B and C are (before rotation):\n";
                        cout << Verts3[0].transpose() << endl;
                        cout << Verts5[0].transpose() << endl;
                        cout << "\nThe two points, B and C are (after rotation):\n";
                        cout << Verts4[0].transpose() << endl;
                        cout << Verts6[0].transpose() << endl;
                        exit(0);
                    }
                }
            }

            //cout << Verts4[0].transpose() << endl;
            //cout << Verts6[0].transpose() << endl;
            if ((Verts6[0](0) > -0.0001 && Verts4[0](0) - Verts6[0](0) > -0.0001) ||
                (Verts6[0](0) < 0.0001 && Verts6[0](0) - Verts4[0](0) > -0.0001))
            {

                EFF = true;
                break;
            }
        }
    }

    return EFF;
};

inline bool Point_2D::If_overlap_with_another_pnt(Point_2D AN)
{
    Vector2d s;
    s = this->Coordinate - AN.Coordinate;
    if (abs(s(0)) < 0.0001 && abs(s(1)) < 0.0001)
    {
        return true;
    }
    else
        return false;
};

inline bool Point_2D::If_lies_on_a_line_seg(const std::vector<Vector2d> Line_seg)
{
    double line_ymax = Line_seg[0](1) > Line_seg[1](1) ? Line_seg[0](1) : Line_seg[1](1);
    double line_ymin = Line_seg[0](1) < Line_seg[1](1) ? Line_seg[0](1) : Line_seg[1](1);

    double line_xmax = Line_seg[0](0) > Line_seg[1](0) ? Line_seg[0](0) : Line_seg[1](0);
    double line_xmin = Line_seg[0](0) < Line_seg[1](0) ? Line_seg[0](0) : Line_seg[1](0);

    // a vertical line, same x _ coordinate
    if (abs(Line_seg[0](0) - Line_seg[1](0)) < 0.0001 && abs(Line_seg[0](1) - Line_seg[1](1)) > 0.0001)
    {
        if (abs(this->Coordinate(0) - Line_seg[0](0)) < 0.0001)
        {
            if (this->Coordinate(1) - line_ymin > -0.0001 && line_ymax - this->Coordinate(1) > -0.0001)
            {
                return true;
            }
            else
                return false;
        }
        else
        {
            return false;
        }
    }

    // a horizontal line
    if (abs(Line_seg[0](0) - Line_seg[1](0)) > 0.0001 && abs(Line_seg[0](1) - Line_seg[1](1)) < 0.0001)
    {
        if (abs(this->Coordinate(1) - Line_seg[0](1)) < 0.0001)
        {
            if (this->Coordinate(0) - line_xmin > -0.0001 && line_xmax - this->Coordinate(0) > -0.0001)
            {
                return true;
            }
            else
                return false;
        }
        else
        {
            return false;
        }
    }

    // a inclined line segment
    Vector2d F = Line_seg[1] - Line_seg[0];
    Vector2d P = this->Coordinate - Line_seg[0];
    std::vector<Vector3d> Verts1(1), Verts2(1), Verts3(1), Verts4(1);
    Verts1[0] << F(0), F(1), 0;
    Verts2[0] << P(0), P(1), 0;

    double angle_line = atan2(F(1), F(0));
    Quaternion_t Q_axis_1;
    Vector3d temp1;
    temp1 << 0, 0, -1;
    DFN::Rotation_verts Re(Verts1, angle_line, Q_axis_1, temp1, Verts3);
    DFN::Rotation_verts Re2(Verts2, angle_line, Q_axis_1, temp1, Verts4);

    if (abs(Verts4[0](1)) > 0.001)
    {
        return false;
    }
    else
    {
        if ((Verts4[0](0) > -0.001 && Verts3[0](0) - Verts4[0](0) > -0.001) ||
            (Verts4[0](0) < 0.001 && Verts4[0](0) - Verts3[0](0) > -0.001))
        {
            return true;
        }
        else
            return false;
    }

    return false;
};

inline Vector2d Point_2D::Perpend_foot_on_a_line_seg(const std::vector<Vector2d> Line_seg)
{
    Vector2d F = Line_seg[1] - Line_seg[0]; // move to origin
    Vector2d P = this->Coordinate - Line_seg[0];
    std::vector<Vector3d> Verts1(1), Verts2(1), Verts3(1), Verts4(1);
    Verts1[0] << F(0), F(1), 0;
    Verts2[0] << P(0), P(1), 0;

    double angle_line = atan2(F(1), F(0));
    Quaternion_t Q_axis_1;
    Vector3d temp1;
    temp1 << 0, 0, -1;
    DFN::Rotation_verts Re(Verts1, angle_line, Q_axis_1, temp1, Verts3); // rotation
    DFN::Rotation_verts Re2(Verts2, angle_line, Q_axis_1, temp1, Verts4);

    std::vector<Vector3d> Verts5(1), Verts6(1);
    Verts5[0] << Verts4[0](0), 0, 0;

    //rotate back
    DFN::Rotation_verts Re3(Verts5, -angle_line, Q_axis_1, temp1, Verts6);

    //move back
    return (Vector2d{Verts6[0](0), Verts6[0](1)} + Line_seg[0]);

};
}; // namespace DFN