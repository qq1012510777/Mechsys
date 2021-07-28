#pragma once
#include "Intersection_between_polygon.h"
#include "Polygon_convex_3D.h"
#include "Project_3D_polygon_to_2D_plane.h"

//Vector6d Model_domain;
///< Top-zmax, bottom-zmin, front-ymin, back-ymax, left-xmin, right-xmax

namespace DFN
{

class Intersection_between_polygon_and_3D_box
{
    std::vector<Vector2d> XY_plane;
    std::vector<Vector2d> XZ_plane;
    std::vector<Vector2d> YZ_plane;

public:
    Intersection_between_polygon_and_3D_box();
    Intersection_between_polygon_and_3D_box(std::vector<Vector3d> &polygon, Vector6d model_domain);
};

inline Intersection_between_polygon_and_3D_box::Intersection_between_polygon_and_3D_box(){
    //
};

inline Intersection_between_polygon_and_3D_box::Intersection_between_polygon_and_3D_box(std::vector<Vector3d> &polygon, Vector6d model_domain)
{
    double Top_zmax = model_domain[0],
           bottom_zmin = model_domain[1],
           front_ymin = model_domain[2],
           back_ymax = model_domain[3],
           left_xmin = model_domain[4],
           right_xmax = model_domain[5];

    XY_plane.resize(4);
    XZ_plane.resize(4);
    YZ_plane.resize(4);

    XY_plane[0] << left_xmin, front_ymin;
    XY_plane[1] << left_xmin, back_ymax;
    XY_plane[2] << right_xmax, back_ymax;
    XY_plane[3] << right_xmax, front_ymin;

    XZ_plane[0] << left_xmin, bottom_zmin;
    XZ_plane[1] << left_xmin, Top_zmax;
    XZ_plane[2] << right_xmax, Top_zmax;
    XZ_plane[3] << right_xmax, bottom_zmin;

    YZ_plane[0] << front_ymin, bottom_zmin;
    YZ_plane[1] << front_ymin, Top_zmax;
    YZ_plane[2] << back_ymax, Top_zmax;
    YZ_plane[3] << back_ymax, bottom_zmin;

    DFN::Polygon_convex_3D poly1{polygon};
    double a = poly1.Plane_parameter[0],
           b = poly1.Plane_parameter[1],
           c = poly1.Plane_parameter[2],
           d = poly1.Plane_parameter[3];

    //check if the normal of polygon is parallel to XY
    Vector3d normal_XY;
    normal_XY << 0, 0, 1;
    if (abs(poly1.Normal_vector.dot(normal_XY)) > 1e-5) // the poly is not perpendicular to XY
    {
        DFN::Project_3D_polygon_to_2D_plane ProjXY;
        std::vector<Vector2d> polygon_XY = ProjXY.Project_to_XY_plane(polygon);

        DFN::Intersection_between_polygon Intersec{polygon_XY, XY_plane};
        if (Intersec.Intersection.size() == 0)
        {
            polygon.clear();
            return;
        }

        polygon_XY = Intersec.Intersection;

        polygon.resize(polygon_XY.size());

        for (size_t i = 0; i < polygon_XY.size(); ++i)
        {
            polygon[i][0] = polygon_XY[i][0];
            polygon[i][1] = polygon_XY[i][1];
            polygon[i][2] = -(a * polygon_XY[i][0] + b * polygon_XY[i][1] + d) / c;
        }
    }

    //check if the normal of polygon is parallel to XZ
    Vector3d normal_XZ;
    normal_XZ << 0, 1, 0;
    if (abs(poly1.Normal_vector.dot(normal_XZ)) > 1e-5) // the poly is not perpendicular to XZ
    {
        DFN::Project_3D_polygon_to_2D_plane ProjXZ;
        std::vector<Vector2d> polygon_XZ = ProjXZ.Project_to_XZ_plane(polygon);

        DFN::Intersection_between_polygon Intersec{polygon_XZ, XZ_plane};
        if (Intersec.Intersection.size() == 0)
        {
            polygon.clear();
            return;
        }

        polygon_XZ = Intersec.Intersection;

        polygon.resize(polygon_XZ.size());

        for (size_t i = 0; i < polygon_XZ.size(); ++i)
        {
            polygon[i][0] = polygon_XZ[i][0];
            polygon[i][1] = -(a * polygon_XZ[i][0] + c * polygon_XZ[i][1] + d) / b;
            polygon[i][2] = polygon_XZ[i][1];
        }
    }

    //check if the normal of polygon is parallel to YZ
    Vector3d normal_YZ;
    normal_YZ << 1, 0, 0;
    if (abs(poly1.Normal_vector.dot(normal_YZ)) > 1e-5) // the poly is not perpendicular to YZ
    {

        DFN::Project_3D_polygon_to_2D_plane ProjYZ;
        std::vector<Vector2d> polygon_YZ = ProjYZ.Project_to_YZ_plane(polygon);

        DFN::Intersection_between_polygon Intersec{polygon_YZ, YZ_plane};
        if (Intersec.Intersection.size() == 0)
        {
            polygon.clear();
            return;
        }

        polygon_YZ = Intersec.Intersection;

        polygon.resize(polygon_YZ.size());

        for (size_t i = 0; i < polygon_YZ.size(); ++i)
        {
            polygon[i][0] = -(b * polygon_YZ[i][0] + c * polygon_YZ[i][1] + d) / a;
            polygon[i][1] = polygon_YZ[i][0];
            polygon[i][2] = polygon_YZ[i][1];
        }
    }
};

}; // namespace DFN