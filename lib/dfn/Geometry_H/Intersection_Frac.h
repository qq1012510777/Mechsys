#pragma once
#include "../DFN_H/Fracture_WL.h"
#include "Convex_3D_polygon_cross_xy_plane.h"
#include "Line_seg_2D.h"
#include "NorVec_plane.h"
#include "Parallel_Inf_Plane.h"
#include "Point_2D.h"
#include "Polygon_convex_3D.h"
#include "Rotation_verts.h"
#include "Vector_2.h"

namespace DFN
{

class Intersection_Frac
{
public:
    bool If_intersect;
    std::vector<Vector3d> Intersection;

public:
    Intersection_Frac(const Polygon_convex_3D F1, const Polygon_convex_3D F2);
};

inline Intersection_Frac::Intersection_Frac(const Polygon_convex_3D F1, const Polygon_convex_3D F2)
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
                std::cout << "Error!!! The Z values of all vertexes of 1st fracture should be the same!\n";
                exit(0);
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
                //exit(0);
            };

            vector<Vector2d> Intersection_1;
            // intersection between Cp.Intersection and the first polygon
            A_r.Intersection_between_a_line_and_a_convex_polygon(Verts_1, Intersection_1);
            //cout << "\ttext 1 \n";
            //cout << "Intersection_1.size(): " << Intersection_1.size() << endl;
            if (Intersection_1.size() == 2)
            {
                //it is a line segment
                //cout << "2" << endl;
                std::vector<Vector3d> PO(2), PT(2);
                PO[0] << Intersection_1[0](0), Intersection_1[0](1), 0;
                PO[1] << Intersection_1[1](0), Intersection_1[1](1), 0;
                PO[0](2) += z_plane;
                PO[1](2) += z_plane;

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
                    PT[1] = PO[1];
                }

                Intersection.resize(2);
                Intersection[0] = PT[0];
                Intersection[1] = PT[1];

                If_intersect = true;
                //cout << "\ttext 2 \n";
                return;
            }
            else if (Intersection_1.size() == 1)
            {

                Vector2d t1, t2;
                t1 << Cp.Intersection[0](0), Cp.Intersection[0](1);
                t2 << Cp.Intersection[1](0), Cp.Intersection[1](1);
                DFN::Line_seg_2D A(t1, t2);

                if (A.If_is_a_point() == true)
                {
                    cout << "Error! found A (Cp.Intersection) is a point\n";
                    exit(0);
                };

                //cout << "\tCp.Intersection:\n";
                //cout << "\t" << Cp.Intersection[0].transpose() << endl;
                //cout << "\t" << Cp.Intersection[1].transpose() << endl;
                DFN::Line_seg_2D NYY{};
                NYY.Extending_this_line_seg(A, Verts_1);

                if (NYY.If_is_a_point() == true)
                {
                    cout << "Error! found NYY (extending Cp.Intersection) is a point\n";
                    exit(0);
                };

                vector<Vector2d> Intersection_A;
                NYY.Intersection_between_a_line_and_a_convex_polygon(Verts_1, Intersection_A);

                //cout << "\tIntersection_A:\n";
                //cout << "\tIntersection_A.size(): " << Intersection_A.size() << endl;
                //cout << "\t" << Intersection_A[0].transpose() << endl;
                //cout << "\t" << Intersection_A[1].transpose() << endl;
                if (Intersection_A.size() == 1)
                {
                    //both Cp.Intersection and Intersection_A touch one of the corner
                    std::vector<Vector3d> PO(1), PT(1);
                    PO[0] << Intersection_A[0](0), Intersection_A[0](1), 0;
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

                    Intersection.resize(2);
                    Intersection[0] = PT[0];
                    Intersection[1] = PT[0];
                    If_intersect = true;
                    return;
                }
                else if (Intersection_A.size() == 2)
                {

                    Vector2d s1, s2;
                    s1 << Intersection_A[0](0), Intersection_A[0](1);
                    s2 << Intersection_A[1](0), Intersection_A[1](1);
                    DFN::Line_seg_2D B(s1, s2);

                    if (B.If_is_a_point() == true)
                    {
                        cout << "Error! found B (Intersection_A between extending line of Cp.Intersection and 1st polygon) is a point\n";
                        exit(0);
                    };

                    vector<Vector2d> Intersection_2;
                    bool dr = A.Intersection_between_two_lines(B, Intersection_2);

                    //cout << "\tIntersection_2.size(): " << Intersection_2.size() << endl;
                    if (Intersection_2.size() == 2)
                    {

                        //it is a line segment
                        std::vector<Vector3d> PO(2), PT(2);
                        PO[0] << Intersection_2[0](0), Intersection_2[0](1), 0;
                        PO[1] << Intersection_2[1](0), Intersection_2[1](1), 0;
                        PO[0](2) += z_plane;
                        PO[1](2) += z_plane;

                        //cout << "\tIntersection_2:\n";
                        //cout << "\t" << Intersection_2[0].transpose() << endl;
                        //cout << "\t" << Intersection_2[1].transpose() << endl;

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
                            PT[1] = PO[1];
                        }

                        Intersection.resize(2);
                        Intersection[0] = PT[0];
                        Intersection[1] = PT[1];
                        If_intersect = true;
                        //cout << "\ttext 2 \n";
                        return;
                    }
                    else if (Intersection_2.size() == 1)
                    {
                        // the Cp.Intersection touches one edge of the 1st polygon from outward,
                        // and so, the Intersection_2 intersect the Cp.Intersection at one point
                        //it is a point
                        std::vector<Vector3d> PO(1), PT(1);
                        PO[0] << Intersection_2[0](0), Intersection_2[0](1), 0;
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

                        Intersection.resize(2);
                        Intersection[0] = PT[0];
                        Intersection[1] = PT[0];
                        If_intersect = true;
                        return;
                    }
                    else
                    {
                        cout << "Intersection between two overlapping lines segments should have one or two ends!!!\n";
                        //cout << "Intersection size: " << Intersection_2.size() << endl;
                        cout << "Line_A=[\n";
                        cout << A.Point[0].transpose() << endl;
                        cout << A.Point[1].transpose() << endl;
                        cout << "];\n";
                        cout << "Line_B=[\n";
                        cout << B.Point[0].transpose() << endl;
                        cout << B.Point[1].transpose() << endl;
                        cout << "];\n";
                        cout << "bool dr: " << dr << endl;

                        cout << "testing if the four points are collinear?\n";
                        bool collinear1, collinear2;
                        DFN::Collinear_2D AS1{A.Point[0], A.Point[1], B.Point[0], collinear1};
                        DFN::Collinear_2D AS2{A.Point[0], A.Point[1], B.Point[1], collinear2};
                        cout << collinear1 << ", " << collinear2 << endl;
                        exit(0);
                    }
                }
                else
                {
                    cout << "Error! extending line of Cp.Intersection  must intersection the 1st polygon!\n";
                    cout << "Because Cp.Intersection intersects the 1st polygon at one point!\n";
                    cout << "Extending line:\n";
                    for (size_t i = 0; i < NYY.Point.size(); ++i)
                        cout << NYY.Point[i].transpose() << endl;
                    cout << "Cp.Intersection (size = " << Cp.Intersection.size() << "):\n";
                    for (size_t i = 0; i < Cp.Intersection.size(); ++i)
                        cout << Cp.Intersection[i].transpose() << endl;

                    cout << "\n1st polygon:\n";
                    for (size_t i = 0; i < Verts_1.size(); ++i)
                        cout << Verts_1[i].transpose() << endl;

                    cout << "\n2nd polygon:\n";
                    for (size_t i = 0; i < Verts_2.size(); ++i)
                        cout << Verts_2[i].transpose() << endl;
                    exit(0);
                }
            }
            else
            {
                //cout << "1" << endl;
                // two potential cases, (1) Cp.Intersection is inside 1st polygon, (2) outside
                Vector2d t1, t2;
                t1 << Cp.Intersection[0](0), Cp.Intersection[0](1);
                t2 << Cp.Intersection[1](0), Cp.Intersection[1](1);
                DFN::Line_seg_2D A(t1, t2);

                if (A.If_is_a_point() == true)
                {
                    cout << "Error! found A (Cp.Intersection) is a point\n";
                    exit(0);
                };

                DFN::Line_seg_2D NYY{};
                NYY.Extending_this_line_seg(A, Verts_1);

                if (NYY.If_is_a_point() == true)
                {
                    cout << "Error! found NYY (extending Cp.Intersection) is a point\n";
                    exit(0);
                };

                vector<Vector2d> Intersection_A;
                NYY.Intersection_between_a_line_and_a_convex_polygon(Verts_1, Intersection_A);

                //cout << "Intersection_A.size(): " << Intersection_A.size() << endl;
                if (Intersection_A.size() == 0)
                {
                    //outside
                    If_intersect = false;
                    return; // no intersection between 1 and 2 polygon
                }
                else if (Intersection_A.size() == 2)
                {
                    Vector2d sq1, sq2;
                    sq1 << Intersection_A[0](0), Intersection_A[0](1);
                    sq2 << Intersection_A[1](0), Intersection_A[1](1);
                    DFN::Line_seg_2D QZ(sq1, sq2);
                    if (QZ.If_is_a_point() == true)
                    {
                        cout << "Error! found QZ (Intersection_A between extending line of Cp.Intersection and 1st polygon) is a point\n";
                        cout << "QZ:\n";
                        cout << QZ.Point[0].transpose() << endl;
                        cout << QZ.Point[1].transpose() << endl;
                        exit(0);
                    };

                    vector<Vector2d> Intersection_B;
                    bool rew = A.Intersection_between_two_lines(QZ, Intersection_B);

                    if (rew == true)
                    {
                        std::vector<Vector3d> PO(2), PT(2);
                        PO[0] << Intersection_B[0](0), Intersection_B[0](1), 0;
                        PO[1] << Intersection_B[1](0), Intersection_B[1](1), 0;
                        PO[0](2) += z_plane;
                        PO[1](2) += z_plane;

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
                            PT[1] = PO[1];
                        }

                        Intersection.resize(2);
                        Intersection[0] = PT[0];
                        Intersection[1] = PT[1];
                        If_intersect = true;
                        //cout << "\ttext 2 \n";
                        return;
                    }
                    else
                    {
                        If_intersect = false;
                        return;
                    }
                }
                else if (Intersection_A.size() == 1)
                {
                    //outside
                    // firstly, Cp.Intersection does not intersect 1st polygon
                    // then the extending of Cp.Intersection intersect 1st polygon at one point
                    // must be non-intersecting
                    If_intersect = false;
                    return; // no intersection between 1 and 2 polygon
                }
            }
        }
        else if (Cp.Relation == "Not_crossing")
        {
            If_intersect = false;
            return;
        }
    }
    else
    {
        cout << "Undefine relation betweem two infinite planes!\n";
        cout << P.Raltion << endl;
        exit(0);
    }
    If_intersect = false;
    return;
}

}; // namespace DFN