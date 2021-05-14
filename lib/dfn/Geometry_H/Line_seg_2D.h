#pragma once
#include "Collinear_2D.h"
#include "Dense"
#include "Interval_1D.h"
#include "Point_2D.h"
#include "Polygon_convex_2D.h"
#include "Rotation_verts.h"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Line_seg_2D
{
public:
    vector<Vector2d> Point;
    Vector2d directional_vector;
    double x_min;
    double x_max;
    double y_min;
    double y_max;

public:
    Line_seg_2D();
    void Re_constructor(const Vector2d P1, const Vector2d P2);
    Line_seg_2D(const Vector2d P1, const Vector2d P2);
    bool Intersection_between_two_lines(const Line_seg_2D Bs, std::vector<Vector2d> &Intersection);
    void Intersection_between_a_line_and_a_convex_polygon(const std::vector<Vector3d> Verts_1, vector<Vector2d> &Intersection);
    void Extream_value();
    bool If_contains_a_polygon(const std::vector<Vector3d> Verts_1);
    void Extending_this_line_seg(const Line_seg_2D Bs, const std::vector<Vector3d> Verts_1);
    bool If_two_lines_overlap(const Line_seg_2D Bs, std::vector<Vector2d> &Intersection);
    bool If_is_a_point();
    vector<Vector2d> Exact_pnts_along_this_line_seg_evenly(const double sub_seg);

private:
    
};

inline Line_seg_2D::Line_seg_2D()
{
    ;
};

inline void Line_seg_2D::Re_constructor(const Vector2d P1, const Vector2d P2)
{
    Point.resize(2);
    Point[0] = P1;
    Point[1] = P2;
    directional_vector = P2 - P1;
    x_min = 0;
    x_max = 0;
    y_min = 0;
    y_max = 0;
};

inline Line_seg_2D::Line_seg_2D(const Vector2d P1, const Vector2d P2)
{
    Point.resize(2);
    Point[0] = P1;
    Point[1] = P2;
    directional_vector = P2 - P1;
    x_min = 0;
    x_max = 0;
    y_min = 0;
    y_max = 0;
};

inline bool Line_seg_2D::Intersection_between_two_lines(const Line_seg_2D Bs, std::vector<Vector2d> &Intersection)
{
    Vector2d T1, T2, U1, U2;
    T1 << this->Point[0](0), this->Point[0](1);
    T2 << this->Point[1](0), this->Point[1](1);
    U1 << Bs.Point[0](0), Bs.Point[0](1);
    U2 << Bs.Point[1](0), Bs.Point[1](1);

    bool itr1, itr2;

    DFN::Collinear_2D KLL1{T1, T2, U1, itr1};
    DFN::Collinear_2D KLL2{T1, T2, U2, itr2};

    if (itr1 == false || itr2 == false)
    {
        segment_type A(point_type(this->Point[0](0), this->Point[0](1)),
                       point_type(this->Point[1](0), this->Point[1](1)));

        segment_type B(point_type(Bs.Point[0](0), Bs.Point[0](1)),
                       point_type(Bs.Point[1](0), Bs.Point[1](1)));

        bool result = boost::geometry::intersects(A, B);
        //cout << "result = " << result << endl;
        if (result == true)
        {
            std::vector<point_type> output;
            boost::geometry::intersection(A, B, output);
            Intersection.resize(output.size());
            for (size_t i = 0; i < Intersection.size(); ++i)
            {
                Intersection[i](0) = output[i].x();
                Intersection[i](1) = output[i].y();
            }
            return true;
        }
        return false;
    }
    else
    {

        // if two lines are overlapping
        std::vector<Vector2d> Intersection_A;
        bool yi = If_two_lines_overlap(Bs, Intersection_A);
        //cout << "yi = " << yi << endl;
        if (yi == false)
        {
            return false;
        }
        else
        {
            Intersection = Intersection_A;
            return true;
        }
    }

    return false;
};

inline bool Line_seg_2D::If_two_lines_overlap(const Line_seg_2D Bs, std::vector<Vector2d> &Intersection)
{
    Vector2d Ay, By, Cy, Dy;
    Ay << this->Point[0](0), this->Point[0](1);
    By << this->Point[1](0), this->Point[1](1);
    Cy << Bs.Point[0](0), Bs.Point[0](1);
    Dy << Bs.Point[1](0), Bs.Point[1](1);

    bool T1, T2;
    DFN::Collinear_2D At{Ay, By, Cy, T1};
    DFN::Collinear_2D Bt{Ay, By, Dy, T2};
    //cout << T1 << T2 << endl;
    if (T1 == true && T2 == true)
    {

        // the two lines must be collinear and overlapping
        Vector2d A1, A2, B1, B2;
        A1 << 0, 0;
        A2 = By - Ay;
        B1 = Cy - Ay;
        B2 = Dy - Ay;

        std::vector<Vector3d> Verts3(3), Verts4(3);
        Verts3[0] << A2(0), A2(1), 0;
        Verts3[1] << B1(0), B1(1), 0;
        Verts3[2] << B2(0), B2(1), 0;

        double angle_tmp = atan2(A2(1), A2(0));
        Quaternion_t Q_axis_1;
        Vector3d temp1;
        temp1 << 0, 0, -1;

        if (abs(Verts3[0](0)) < 0.0001 && abs(Verts3[1](0)) < 0.0001 && abs(Verts3[2](0)) < 0.0001)
        {
            angle_tmp = 90.0 * M_PI / 180.0;
        }

        if (abs(A2(1)) > 0.0001 || abs(B1(1)) > 0.0001 || abs(B2(1)) > 0.0001)
            DFN::Rotation_verts Re(Verts3, angle_tmp, Q_axis_1, temp1, Verts4);
        else
            Verts4 = Verts3;

        if (abs(Verts4[0](1)) > 0.0001 ||
            abs(Verts4[1](1)) > 0.0001 ||
            abs(Verts4[2](1)) > 0.0001)
        {
            cout << "Line_seg_2D::If_two_lines_overlap, Rotating to x_axis failed!\n";
            cout << "three y values: " << Verts4[0](1) << ", " << Verts4[1](1) << ", " << Verts4[2](1) << endl
                 << endl;
            cout << "the first line:\n";
            cout << this->Point[0].transpose() << endl;
            cout << this->Point[1].transpose() << endl;
            cout << "the second line:\n";
            cout << Bs.Point[0].transpose() << endl;
            cout << Bs.Point[1].transpose() << endl;

            cout << "\nafter move, first line:\n";
            cout << A1.transpose() << endl;
            cout << A2.transpose() << endl;
            cout << "\nafter move, second line:\n";
            cout << B1.transpose() << endl;
            cout << B2.transpose() << endl;

            cout << "\nafter rotation, first line:\n";
            cout << A1.transpose() << endl;
            cout << Verts4[0].transpose() << endl;
            cout << "\nafter rotation, second line:\n";
            cout << Verts4[1].transpose() << endl;
            cout << Verts4[2].transpose() << endl;

            cout << "angle of 1st line: " << angle_tmp * 180. / M_PI << endl;
            exit(0);
        }

        //now they are rotated to be adhering to x-axis
        Vector2d ASS, BSS;
        ASS << 0, Verts4[0](0);
        BSS << Verts4[1](0), Verts4[2](0);

        DFN::Interval_1D interval{ASS, BSS};
        /*
        cout << "ASS: " << ASS.transpose() << endl;
        cout << "BSS: " << BSS.transpose() << endl;
        cout << "\n";
        cout << "ASS: " << interval.x1_min << ", " << interval.x1_max << endl;
        cout << "BSS: " << interval.x2_min << ", " << interval.x2_max << endl;
        cout << "Intersection 1D: " << interval.Intersection.transpose() << endl;
        */
        if (interval.if_intersect == true)
        {

            std::vector<Vector3d> Verts5(2), Verts6(2);

            for (size_t i = 0; i < Verts5.size(); ++i)
            {
                Verts5[i](0) = interval.Intersection(i);
                Verts5[i](1) = 0;
                Verts5[i](2) = 0;
            }

            //rotation back
            if (abs(angle_tmp * 180. / M_PI) > 0.0001)
                DFN::Rotation_verts Re(Verts5, -angle_tmp, Q_axis_1, temp1, Verts6);
            else
                Verts6 = Verts5;
            //cout << "\tVerts6:\n";
            //cout << "\t" << Verts6[0].transpose() << endl;
            //cout << "\t" << Verts6[1].transpose() << endl;
            Intersection.resize(Verts6.size());
            for (size_t i = 0; i < Intersection.size(); ++i)
            {
                Intersection[i] << Verts6[i](0) + Ay(0), Verts6[i](1) + Ay(1);
            }

            return true;
        }
        else
            return false;
    }
    else
    {
        return false;
    }
};

inline void Line_seg_2D::Intersection_between_a_line_and_a_convex_polygon(const std::vector<Vector3d> Verts_1, vector<Vector2d> &Intersection)
{
    for (size_t i = 0; i < Verts_1.size(); ++i)
    {
        Vector2d t1, t2;
        t1 << Verts_1[i](0), Verts_1[i](1);
        t2 << Verts_1[(i + 1) % Verts_1.size()](0), Verts_1[(i + 1) % Verts_1.size()](1);

        DFN::Line_seg_2D TMP_z{t1, t2};
        std::vector<Vector2d> IN_A;
        bool tu = Intersection_between_two_lines(TMP_z, IN_A);

        if (tu == true && IN_A.size() == 1)
        {
            Intersection.push_back(IN_A[0]);
        }
        else if (tu == true && IN_A.size() == 2)
        {
            Intersection.resize(2);
            Intersection[0] = IN_A[0];
            Intersection[1] = IN_A[1];
            return;
        }
    }

    if (Intersection.size() >= 2 && Intersection.size() <= 4)
    {
        //a special case exists that the line intersects one or two corners of the polygon, which yields two to four intersections
        vector<Vector2d> Intersection_A;
        Intersection_A.push_back(Intersection[0]);
        for (size_t i = 1; i < Intersection.size(); ++i)
        {
            bool aka = false;
            for (size_t j = i - 1;; j--)
            {
                if (abs(Intersection[i](0) - Intersection[j](0)) < 0.0001 ||
                    abs(Intersection[i](1) - Intersection[j](1)) < 0.0001)
                {
                    aka = true;
                    break;
                }
                if (j == 0)
                {
                    break;
                }
            }

            if (aka == false)
            {
                Intersection_A.push_back(Intersection[i]);
            }
        }

        if (Intersection_A.size() > 2)
        {
            cout << "Line_seg_2D::Intersection_between_a_line_and_a_convex_polygon, a line cannot intersect a convex polygon more than two times!\n";
            exit(0);
        }

        Intersection.resize(Intersection_A.size());
        for (size_t i = 0; i < Intersection.size(); ++i)
        {
            Intersection[i] = Intersection_A[i];
        }
    }
    else if (Intersection.size() > 4)
    {
        cout << "\tA line segment cannot intersect a convex polygon more than two times!\n";
        cout << "\tIntersection is:\n";
        for (size_t i = 0; i < Intersection.size(); ++i)
            cout << "\t" << Intersection[i].transpose() << endl;
        exit(0);
    }
    else if (Intersection.size() <= 1)
    {
        ;
    }
};

inline void Line_seg_2D::Extream_value()
{
    if (Point[0](0) > Point[1](0))
    {
        x_min = Point[1](0);
        x_max = Point[0](0);
    }
    else
    {
        x_min = Point[0](0);
        x_max = Point[1](0);
    }

    if (Point[0](1) > Point[1](1))
    {
        y_min = Point[1](1);
        y_max = Point[0](1);
    }
    else
    {
        y_min = Point[0](1);
        y_max = Point[1](1);
    }
};

inline bool Line_seg_2D::If_contains_a_polygon(const std::vector<Vector3d> Verts_1)
{
    Extream_value();
    DFN::Polygon_convex_2D polygon1(Verts_1);
    polygon1.Extream();
    if ((x_min < polygon1.x_min &&
         x_max > polygon1.x_max) ||
        (y_min < polygon1.y_min &&
         y_max > polygon1.y_max))
    {
        return true;
    }
    else
        return false;
};

inline void Line_seg_2D::Extending_this_line_seg(const Line_seg_2D Bs, const std::vector<Vector3d> Verts_1)
{
    Vector2d A, B;

    if (abs(Bs.directional_vector(0)) > 0.0001 && abs(Bs.directional_vector(1)) > 0.0001)
    {
        for (size_t i = 0;; i++)
        {
            double x1 = (i + 1) * (-100.);
            double y1 = (x1 - Bs.Point[0](0)) / Bs.directional_vector(0) * Bs.directional_vector(1) + Bs.Point[0](1);

            double x2 = (i + 1) * (100.);
            double y2 = (x2 - Bs.Point[0](0)) / Bs.directional_vector(0) * Bs.directional_vector(1) + Bs.Point[0](1);

            A << x1, y1;
            B << x2, y2;

            DFN::Line_seg_2D New(A, B);
            if (New.If_contains_a_polygon(Verts_1) == true)
            {
                break;
            }
        }
    }
    else if (abs(Bs.directional_vector(0)) < 0.0001 && abs(Bs.directional_vector(1)) > 0.0001)
    {
        for (size_t i = 0;; i++)
        {
            double x1 = Bs.Point[1](0), x2 = Bs.Point[1](0);
            double y1 = -100. * (i + 1);
            double y2 = 100. * (i + 1);

            A << x1, y1;
            B << x2, y2;

            DFN::Line_seg_2D New(A, B);
            if (New.If_contains_a_polygon(Verts_1) == true)
            {
                break;
            }
        }
    }
    else if (abs(Bs.directional_vector(0)) > 0.0001 && abs(Bs.directional_vector(1)) < 0.0001)
    {
        for (size_t i = 0;; i++)
        {

            double y1 = Bs.Point[1](1), y2 = Bs.Point[1](1);
            double x1 = -100. * (i + 1);
            double x2 = 100. * (i + 1);

            A << x1, y1;
            B << x2, y2;

            DFN::Line_seg_2D New(A, B);
            if (New.If_contains_a_polygon(Verts_1) == true)
            {
                break;
            }
        }
    }

    Point.resize(2);
    Point[0] = A;
    Point[1] = B;
    directional_vector = B - A;
};

inline bool Line_seg_2D::If_is_a_point()
{
    Vector2d su;
    su = this->Point[0] - this->Point[1];

    if (abs(su(0)) < 0.0001 && abs(su(1)) < 0.0001)
    {
        return true;
    }
    else
        return false;
};

inline vector<Vector2d> Line_seg_2D::Exact_pnts_along_this_line_seg_evenly(const double sub_seg)
{
    vector<Vector2d> fr;

    if (this->If_is_a_point() == true)
    {
        return fr;
    }

    double length = (this->Point[1] - this->Point[0]).norm();
    size_t numofsubseg = round((length / sub_seg), 0);

    if (numofsubseg == 0)
        return fr;

    double a_sub_seg = length / numofsubseg;
    size_t pntnum = numofsubseg - 1;
    fr.resize(pntnum);

    if (pntnum == 0)
    {
        return fr;
    }

    Vector2d A;
    A << this->Point[1] - this->Point[0];

    double angle = atan2(A(1), A(0));
    Vector3d z_axis{0, 0, -1};

    Quaternion_t Q_axis;
    std::vector<Vector3d> Verts1(1), Verts2(1);
    Verts1[0] << A(0), A(1), 0;

    DFN::Rotation_verts Re(Verts1, angle, Q_axis, z_axis, Verts2);

    std::vector<Vector3d> Verts3(pntnum), Verts4(pntnum);
    for (size_t i = 0; i < pntnum; ++i)
    {
        Verts3[i] << (i + 1) * a_sub_seg, 0, 0;
    }

    DFN::Rotation_verts Rew(Verts3, -angle, Q_axis, z_axis, Verts4);

    for (size_t i = 0; i < pntnum; ++i)
    {
        fr[i] = Vector2d{Verts4[i](0), Verts4[i](1)} + this->Point[0];
    }

    return fr;
};

}; // namespace DFN