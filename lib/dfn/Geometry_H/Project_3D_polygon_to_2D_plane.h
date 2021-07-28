#pragma once
#include "Polygon_convex_2D.h"

namespace DFN
{

class Project_3D_polygon_to_2D_plane
{
public:
    Project_3D_polygon_to_2D_plane();
    std::vector<Vector2d> Project_to_XY_plane(std::vector<Vector3d> polygon);
    std::vector<Vector2d> Project_to_XZ_plane(std::vector<Vector3d> polygon);
    std::vector<Vector2d> Project_to_YZ_plane(std::vector<Vector3d> polygon);
};

inline Project_3D_polygon_to_2D_plane::Project_3D_polygon_to_2D_plane(){
    //
};

inline std::vector<Vector2d> Project_3D_polygon_to_2D_plane::Project_to_XY_plane(std::vector<Vector3d> polygon)
{
    std::vector<Vector2d> XY(polygon.size());
    for (size_t i = 0; i < polygon.size(); ++i)
        XY[i] << polygon[i][0], polygon[i][1];

    return XY;
};

inline std::vector<Vector2d> Project_3D_polygon_to_2D_plane::Project_to_XZ_plane(std::vector<Vector3d> polygon)
{
    std::vector<Vector2d> XZ(polygon.size());
    for (size_t i = 0; i < polygon.size(); ++i)
        XZ[i] << polygon[i][0], polygon[i][2];

    return XZ;
};

inline std::vector<Vector2d> Project_3D_polygon_to_2D_plane::Project_to_YZ_plane(std::vector<Vector3d> polygon)
{
    std::vector<Vector2d> YZ(polygon.size());
    for (size_t i = 0; i < polygon.size(); ++i)
        YZ[i] << polygon[i][1], polygon[i][2];

    return YZ;
};

}; // namespace DFN