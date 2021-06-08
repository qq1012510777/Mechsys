#pragma once
#include "../Math_WL_H/Math_WL.h"
#include "Dense"
#include "Line_seg_3D.h"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Convex_3D_polygon_cross_xy_plane
{
public:
    string Relation;
    std::vector<Vector3d> Intersection;

public:
    Convex_3D_polygon_cross_xy_plane(const std::vector<Vector3d> Verts_2);
};

inline Convex_3D_polygon_cross_xy_plane::Convex_3D_polygon_cross_xy_plane(const std::vector<Vector3d> Verts_2)
{
    size_t Idx_large = 0;
    size_t Idx_small = 0;
    size_t Idx_equal = 0;

    vector<Vector3d> touching_pnt;
    /*
    cout << "std::vector<Vector3d> Verts_2:\n";
    for (size_t i = 0; i < Verts_2.size(); ++i)
        cout << Verts_2[i].transpose() << endl;
    */
    for (size_t i = 0; i < Verts_2.size(); ++i)
    {
        if (round(abs(Verts_2[i](2)), 4) < 0.0001)
        {
            Idx_equal++;
            touching_pnt.push_back(Verts_2[i]);
            //cout << "push: " << Verts_2[i].transpose() << endl;
            goto Y1600;
        }

        if (Verts_2[i](2) > 0.0001)
        {
            Idx_large++;
            goto Y1600;
        }

        if (Verts_2[i](2) < -0.0001)
        {
            Idx_small++;
            goto Y1600;
        }
    Y1600:;
        ;
    }

    //special case 1: polygon intersects xy_plane at a point
    if (Idx_equal == 1 && (Idx_large == 0 || Idx_small == 0))
    {
        Relation = "Intersecting_at_a_point";
        Intersection = touching_pnt;
        return;
    }

    //special case 2: an edge of the polygon lies on the xy_plane
    if (Idx_equal == 2)
    {
        Relation = "An_edge_on_xy_plane";
        Intersection = touching_pnt;
        return;
    }

    if (Idx_equal > 2)
    {
        throw Error_throw_ignore("more than one edge lies on the xy_plane!\nIn class 'Convex_3D_polygon_cross_xy_plane', function 'Convex_3D_polygon_cross_xy_plane'\n");
    }

    if (Idx_large >= 1 && Idx_small >= 1)
    {
        Relation = "Crossing";
        for (size_t i = 0; i < Verts_2.size(); ++i)
        {
            DFN::Line_seg_3D line_seg_3D(Verts_2[i], Verts_2[(i + 1) % Verts_2.size()]);
            Vector3d A;
            bool ty = line_seg_3D.Line_crosses_a_horizontal_plane(A);
            if (ty == true)
            {
                Intersection.push_back(A);
            }
        }

        if (Intersection.size() != 2)
        {
            cout << "Error! in Convex_3D_polygon_cross_xy_plane::Convex_3D_polygon_cross_xy_plane!\n";
            cout << "The 2nd polygon cross the horizontal plabe, but " << Intersection.size() << " ends of intersection were found!\n";
            cout << "\n2nd polygon:\n";
            for (size_t i = 0; i < Verts_2.size(); ++i)
                cout << Verts_2[i].transpose() << endl;

            cout << "\nIntersection:\n";
            for (size_t i = 0; i < Intersection.size(); ++i)
                cout << Intersection[i].transpose() << endl;
            
            throw Error_throw_ignore("Error! in Convex_3D_polygon_cross_xy_plane::Convex_3D_polygon_cross_xy_plane!\n");
        }
        return;
    }

    Relation = "Not_crossing";
};

}; // namespace DFN