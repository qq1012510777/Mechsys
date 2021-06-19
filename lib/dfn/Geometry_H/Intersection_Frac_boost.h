#pragma once
#include "../DFN_H/Fracture_WL.h"
#include "Convex_3D_polygon_cross_xy_plane.h"
#include "Line_seg_2D.h"
#include "NorVec_plane.h"
#include "Parallel_Inf_Plane.h"
#include "Point_2D.h"
#include "Point_3D.h"
#include "Polygon_convex_3D.h"
#include "Rotation_verts.h"
#include "Vector_2.h"

//-- boost
#define BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING

#include "boost/assign/std/vector.hpp"
#include "boost/geometry.hpp"
#include "boost/geometry/algorithms/area.hpp"
#include "boost/geometry/algorithms/assign.hpp"
#include "boost/geometry/algorithms/intersection.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/geometries/polygon.hpp"
#include "boost/geometry/io/dsv/write.hpp"

typedef boost::geometry::model::linestring<point_type> Linestring;
typedef boost::geometry::model::multi_point<point_type> BMultiPoint;
//----------------------------------

namespace DFN
{

class Intersection_Frac_boost
{
public:
    bool If_intersect;
    std::vector<Vector3d> Intersection;

public:
    Intersection_Frac_boost(const Polygon_convex_3D F1, const Polygon_convex_3D F2);
};

inline Intersection_Frac_boost::Intersection_Frac_boost(const Polygon_convex_3D F1, const Polygon_convex_3D F2)
{

    /* std::cout<<"debug\n"; */

    DFN::Parallel_Inf_Plane P(F1, F2);
    //cout << P.Raltion << endl;

    if (P.Raltion == "Parallel")
    {
        If_intersect = false;
        return;
    }
    else if (P.Raltion == "Overlapping")
    {
        If_intersect = false;
        // I have to consider if two fractures are overlapping, how to address!
        // but this happens very very very rarely!!!
        return;
    }
    else if (P.Raltion == "Intersecting")
    {
        vector<Vector3d> f1, f2;
        f1 = F1.Corners;
        f2 = F2.Corners;

        std::vector<Vector3d> Verts_1;
        std::vector<Vector3d> Verts_2;

        Vector3d temp1;
        DFN::Vector_2 v(F1.Normal_vector, temp1);

        double R_angle_temp1 = 0;
        double x_temp = F1.Beta;
        R_angle_temp1 = x_temp * M_PI / 180;
        Quaternion_t Q_axis_1;

        if (F1.Beta > 0.0001)
        {
            //cout << "1;\n";
            Verts_1.resize(f1.size());
            Verts_2.resize(f2.size());
            DFN::Rotation_verts R1(f1, R_angle_temp1, Q_axis_1, temp1, Verts_1);
            DFN::Rotation_verts R2(f2, R_angle_temp1, Q_axis_1, temp1, Verts_2);
        }
        else
        {
            Verts_1 = f1;
            Verts_2 = f2;
        }

        //-------a piece of debug code----
        for (size_t i = 0; i < Verts_1.size(); i++)
        {
            if (abs(Verts_1[i](2) - Verts_1[(i + 1) % Verts_1.size()](2)) > 0.0001)
            {
                DFN::Polygon_convex_3D ty1(Verts_1);
                DFN::Polygon_convex_3D ty2(Verts_2);

                cout << "\npolygon 1\n";
                ty1.Print_Polygon();
                cout << "\npolygon 2\n";
                ty2.Print_Polygon();
                cout << "Normal vector:\n";
                cout << F1.Normal_vector.transpose() << endl;
                cout << "Rotation vector:\n";
                cout << temp1.transpose() << endl;
                //cout << F2.Normal_vector.transpose() << endl;

                cout << "calculated normal vector:\n";
                DFN::NorVec_plane CYT{f1};
                cout << CYT.Normal_vector.transpose() << endl;

                cout << "if parallel: \n";
                cout << F1.Normal_vector(0) / CYT.Normal_vector(0) << ", ";
                cout << F1.Normal_vector(1) / CYT.Normal_vector(1) << ", ";
                cout << F1.Normal_vector(2) / CYT.Normal_vector(2) << "\n\n";

                cout << "if use calculated normal:\n";
                DFN::Vector_2 vs(CYT.Normal_vector, temp1);
                DFN::Rotation_verts R1(f1, R_angle_temp1, Q_axis_1, temp1, Verts_1);
                for (size_t j = 0; j < Verts_1.size(); ++j)
                    cout << Verts_1[j].transpose() << endl;
                throw Error_throw_ignore("Error!!! The Z values of all vertexes of 1st fracture should be the same!\n");
            }
        }
        //--------------------------------

        // move Verts_1 to x-y plane
        for (size_t i = 0; i < Verts_2.size(); i++)
            Verts_2[i](2) -= Verts_1[1](2);

        double z_plane = Verts_1[1](2);
        for (size_t i = 0; i < Verts_1.size(); i++)
            Verts_1[i](2) = 0;

        DFN::Convex_3D_polygon_cross_xy_plane Cp(Verts_2);
        //cout << "\t" << Cp.Relation << endl;
        if (Cp.Relation != "Not_crossing")
        {
            /*
            for (size_t i = 0; i < Cp.Intersection.size(); ++i)
                cout << "\t" << Cp.Intersection[i].transpose() << endl;
            */
            //return false;
        }

        if (Cp.Relation == "Intersecting_at_a_point")
        {

        KF150:;
            Cp.Relation = "Intersecting_at_a_point";
            Vector2d As;
            As << Cp.Intersection[0](0), Cp.Intersection[0](1);
            DFN::Point_2D PNT(As);
            bool tu = PNT.If_lies_within_a_polygon_2D(Verts_1);
            bool ts = PNT.If_lies_on_the_bounds_of_polygon(Verts_1);

            if (tu == true || ts == true)
            {
                std::vector<Vector3d> PO(1), PT(1);
                PO[0] = Cp.Intersection[0];
                PO[0](2) += z_plane;

                Quaternion_t Q_axis_2;
                Vector3d temp2;
                DFN::Vector_2 v2(F1.Normal_vector, temp2);
                if (F1.Beta > 0.0001)
                {
                    DFN::Rotation_verts Rs(PO, -R_angle_temp1, Q_axis_2, temp2, PT);
                }
                else
                {
                    PT[0] = PO[0];
                }

                //cout << "\tPT: " << PT[0].transpose() << endl;
                If_intersect = true;
                Intersection.resize(2);
                Intersection[0] = PT[0];
                Intersection[1] = PT[0];
                return;
            }
        }
        else if (Cp.Relation == "An_edge_on_xy_plane" || Cp.Relation == "Crossing")
        {
            /*
            Cp.Intersection[0](0), Cp.Intersection[0](1);
            Cp.Intersection[1](0), Cp.Intersection[1](1); // they are the trace of polygon 2 on horizontal plane
            */
            Vector2d t1, t2;
            t1 << Cp.Intersection[0](0), Cp.Intersection[0](1);
            t2 << Cp.Intersection[1](0), Cp.Intersection[1](1);

            DFN::Line_seg_2D A_r(t1, t2);
            if (A_r.If_is_a_point() == true)
            {
                /*
                cout << "Error! found Cp.Relation is (An_edge_on_xy_plane) or (Crossing), but actually Cp.Intersection is a point\n";
                cout << "Cp.Relation: " << Cp.Relation << endl;
                cout << "Cp.Intersection:\n";
                cout << Cp.Intersection[0].transpose() << endl;
                cout << Cp.Intersection[1].transpose() << endl;
                cout << "\nAftrt rotating to horizontal plane, 2nd polygon becomes:\n";
                for (size_t k = 0; k < Verts_2.size(); ++k)
                    cout << Verts_2[k].transpose() << endl;
                */
                goto KF150;
               
            };

            std::vector<std::vector<double>> UO(2);
            Linestring TraceOf_2_On_horizon{{Cp.Intersection[0](0), Cp.Intersection[0](1)}, {Cp.Intersection[1](0), Cp.Intersection[1](1)}};

            std::vector<Vector2d> POLYGON_1(Verts_1.size());
            for (size_t i = 0; i < POLYGON_1.size(); ++i)
            {
                POLYGON_1[i] << Verts_1[i](0), Verts_1[i](1);
            }

            std::vector<point_type> points(POLYGON_1.size() + 1);

            for (size_t i = 0; i < POLYGON_1.size(); ++i)
                points[i] = point_type(POLYGON_1[i](0), POLYGON_1[i](1));

            points[points.size() - 1] = point_type(POLYGON_1[0](0), POLYGON_1[0](1));

            polygon_type poly_2;
            boost::geometry::assign_points(poly_2, points);

            DFN::Point_2D ENDA{Vector2d{Cp.Intersection[0](0), Cp.Intersection[0](1)}};
            DFN::Point_2D ENDB{Vector2d{Cp.Intersection[1](0), Cp.Intersection[1](1)}};

            bool A_in = ENDA.If_lies_within_a_polygon_boost(POLYGON_1);
            bool A_on_edge = ENDA.If_lies_on_the_edges_of_a_polygon_boost(POLYGON_1);
            bool B_in = ENDB.If_lies_within_a_polygon_boost(POLYGON_1);
            bool B_on_edge = ENDB.If_lies_on_the_edges_of_a_polygon_boost(POLYGON_1);

            std::vector<Vector3d> Intersection_Two_Polygon;

            if ((A_in == true && B_in == true) ||
                (A_on_edge == true && B_on_edge == true) ||
                (A_in == true && B_on_edge == true) ||
                (A_on_edge == true && B_in == true))
            {
                Intersection_Two_Polygon.resize(2);
                Intersection_Two_Polygon[0] << Cp.Intersection[0](0), Cp.Intersection[0](1), 0;
                Intersection_Two_Polygon[1] << Cp.Intersection[1](0), Cp.Intersection[1](1), 0;
            }
            else if (
                ((A_in == false && A_on_edge == false) && (B_in == true)) ||
                ((A_in == true) && (B_in == false && B_on_edge == false)))
            {
                Intersection_Two_Polygon.resize(2);
                if (B_in == true)
                {
                    Intersection_Two_Polygon[0] << Cp.Intersection[1](0), Cp.Intersection[1](1), 0;
                }

                if (A_in == true)
                {
                    Intersection_Two_Polygon[0] << Cp.Intersection[0](0), Cp.Intersection[0](1), 0;
                }

                std::vector<point_type> output;

                boost::geometry::intersection(TraceOf_2_On_horizon, poly_2, output);
                Intersection_Two_Polygon[1] << output[0].x(), output[0].y(), 0;
            }
            else if (((A_in == false && A_on_edge == false) && (B_on_edge == true)) ||
                     ((A_on_edge == true) && (B_in == false && B_on_edge == false)))
            {
                Intersection_Two_Polygon.resize(1);
                if (B_on_edge == true)
                {
                    Intersection_Two_Polygon[0] << Cp.Intersection[1](0), Cp.Intersection[1](1), 0;
                }

                if (A_on_edge == true)
                {
                    Intersection_Two_Polygon[0] << Cp.Intersection[0](0), Cp.Intersection[0](1), 0;
                }

                std::vector<point_type> output;

                boost::geometry::intersection(TraceOf_2_On_horizon, poly_2, output);

                for (size_t i = 0; i < output.size(); ++i)
                {
                    Vector3d aus;
                    aus << output[i].x(), output[i].y(), 0;

                    if ((aus - Intersection_Two_Polygon[0]).norm() > 0.01)
                    {
                        Intersection_Two_Polygon.push_back(aus);
                    }
                }

                if (Intersection_Two_Polygon.size() > 2)
                {
                    std::vector<DFN::Point_3D> FTK(Intersection_Two_Polygon.size());
                    for (size_t i = 0; i < Intersection_Two_Polygon.size(); ++i)
                        FTK[i].Re_constructor(Intersection_Two_Polygon[i]);
                    
                    DFN::Point_3D TY;

                    TY.Remove_the_same_pnt(FTK);
                   
                    //--------------------------
                    if (FTK.size() > 2)
                    {
                        cout << "Error! the size of Intersection_Two_Polygon cannot more than two!\n";
                        for (size_t i = 0; i < Intersection_Two_Polygon.size(); ++i)
                            cout << Intersection_Two_Polygon[i].transpose() << endl;

                        throw Error_throw_ignore("Error! the size of Intersection_Two_Polygon cannot more than two!\n");
                    }
                    //----------------------------

                    Intersection_Two_Polygon.resize(FTK.size());
                    for (size_t i = 0; i < FTK.size(); ++i)
                        Intersection_Two_Polygon[i] = FTK[i].Coordinate;  
                }
            }
            else // two ends are outside the polygon2
            {
                std::vector<point_type> output;

                boost::geometry::intersection(TraceOf_2_On_horizon, poly_2, output);

                if (output.size() == 0)
                {
                    If_intersect = false;
                    return;
                }

                std::vector<Vector2d> TMP_INTERSECT(output.size());
                for (size_t i = 0; i < output.size(); ++i)
                    TMP_INTERSECT[i] << output[i].x(), output[i].y();

                if (TMP_INTERSECT.size() == 1)
                {
                    Intersection_Two_Polygon.resize(1);
                    Intersection_Two_Polygon[0] << output[0].x(), output[0].y(), 0;
                }
                else if (TMP_INTERSECT.size() == 2)
                {
                    if ((TMP_INTERSECT[0] - TMP_INTERSECT[1]).norm() > 0.01)
                    {
                        Intersection_Two_Polygon.resize(2);
                        Intersection_Two_Polygon[0] << output[0].x(), output[0].y(), 0;
                        Intersection_Two_Polygon[1] << output[1].x(), output[1].y(), 0;
                    }
                    else
                    {
                        Intersection_Two_Polygon.resize(1);
                        Intersection_Two_Polygon[0] << output[0].x(), output[0].y(), 0;
                    }
                }
                else if (TMP_INTERSECT.size() > 2)
                {
                    cout << "Error! A line segment intersect a convex polygon more than 2 times!\nThe line seg is:\n";
                    for (size_t i = 0; i < output.size(); ++i)
                        cout << TMP_INTERSECT[i].transpose() << endl;
                    cout << "The polygon is:\n";
                    for (size_t i = 0; i < POLYGON_1.size(); ++i)
                        cout << POLYGON_1[i].transpose() << endl;
                    throw Error_throw_ignore("Error! A line segment intersect a convex polygon more than 2 times!\n");
                }
                else
                {
                    If_intersect = false;
                    return;
                }
            }
            // now we found Intersection_Two_Polygon , rotate back

            for (size_t i = 0; i < Intersection_Two_Polygon.size(); ++i)
            {
                Intersection_Two_Polygon[i](2) += z_plane;
            }
            std::vector<Vector3d> PT(Intersection_Two_Polygon.size());

            Quaternion_t Q_axis_2;
            Vector3d temp2;
            DFN::Vector_2 v2(F1.Normal_vector, temp2);
            if (F1.Beta > 0.0001)
            {
                DFN::Rotation_verts Rs(Intersection_Two_Polygon, -R_angle_temp1, Q_axis_2, temp2, PT);
            }
            else
            {
                for (size_t i = 0; i < Intersection_Two_Polygon.size(); ++i)
                    PT[i] = Intersection_Two_Polygon[i];
            }

            //cout << "\tPT: " << PT[0].transpose() << endl;
            If_intersect = true;
            Intersection.resize(2);
            if (PT.size() == 2)
            {
                Intersection[0] = PT[0];
                Intersection[1] = PT[1];
            }
            else if (PT.size() == 1)
            {
                Intersection[0] = PT[0];
                Intersection[1] = PT[0];
            }
            return;
        }
        else if (Cp.Relation == "Not_crossing")
        {
            If_intersect = false;
            return;
        }
    }
    else
    {
        cout << P.Raltion << endl;
        throw Error_throw_ignore("Undefine relation betweem two infinite planes!\n");
    }
    If_intersect = false;
    return;
}

}; // namespace DFN