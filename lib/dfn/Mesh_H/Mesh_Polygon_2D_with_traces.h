#pragma once
#include "../Geometry_H/Splitting_Polygon_convex_2D_with_traces.h"
#include "../Math_WL_H/Math_WL.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Fade_2D
#include "Fade_2D.h"
/*_______________________________
here, we added an open source triangulation library, Fade_2D
https://www.geom.at/fade2d/html/
Dr. Bernhard Kornberger
_______________________________*/

using namespace GEOM_FADE2D;
using namespace std;

typedef Eigen::Matrix<int, 1, 6> RowVector6i;

typedef struct If_special_pnt
{
    bool If_frac_bound = false;
    bool If_trace = false;

    /*
    bool If_model_top = false;
    bool If_model_bottom = false;
    bool If_model_front = false;
    bool If_model_back = false;
    bool If_model_left = false;
    bool If_model_right = false;
    */
    //cannot identify if this point is model boundary point or not at this moment, because it is in 2D space
    //Top-zmax, bottom-zmin, front-ymin, back-ymax, left-xmin, right-xmax

} spe_pnt;

namespace DFN
{
//---------------------------------
class YourPeelPredicate_s : public UserPredicateT
{
public:
    bool operator()(const Triangle2 *pT)
    {
        // Return true when the triangle has an interior angle below 1.0
        for (int i = 0; i < 3; ++i)
        {
            if (pT->getInteriorAngle2D(i) < 1.0)
                return true;
        }
        return false;
    }
};
//---------------------------------
class Mesh_Polygon_2D_with_traces
{
public:
    std::vector<RowVector6i> JM;
    std::vector<Eigen::Vector2d> JXY;
    std::vector<spe_pnt> Pnt_attribute;
    std::vector<size_t> Trace_Tag;

public:
    Mesh_Polygon_2D_with_traces(const std::vector<Vector2d> Pnt_sets,
                                const std::vector<pair<size_t, size_t>> Topo_struct,
                                const double min_angle,
                                const double min_edge_tri,
                                const double max_edge_tri,
                                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared,
                                size_t &NO_Nodes_p);

    void Find_neighbor_ele_and_shared_edge(const std::vector<RowVector6i> JM,
                                           std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared);

    void Insert_central_pnt(std::vector<Eigen::Vector2d> &JXY,
                            std::vector<RowVector6i> &JM,
                            const std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> neigh_shared);

    void Matlab_plot(string FileKey, const DFN::Splitting_Polygon_convex_2D_with_traces Poly_t);

    void Identifying_Frac_bound(const DFN::Polygon_convex_2D polygon, const vector<DFN::Line_seg_2D> Traces);
};

inline Mesh_Polygon_2D_with_traces::Mesh_Polygon_2D_with_traces(const std::vector<Vector2d> Pnt_sets,
                                                                const std::vector<pair<size_t, size_t>> Topo_struct,
                                                                const double min_angle,
                                                                const double min_edge_tri,
                                                                const double max_edge_tri,
                                                                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared,
                                                                size_t &NO_Nodes_p)
{
    std::vector<Point2> vInputPoints;
    vInputPoints.resize(Pnt_sets.size());

    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        vInputPoints[i] = Point2(Pnt_sets[i](0), Pnt_sets[i](1));
    }

    Fade_2D dt;
    dt.insert(vInputPoints);
    std::vector<Segment2> vSegments;
    vSegments.resize(Topo_struct.size());
    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        vSegments[i] = Segment2(vInputPoints[Topo_struct[i].first], vInputPoints[Topo_struct[i].second]);
    }

    dt.createConstraint(vSegments, CIS_CONSTRAINED_DELAUNAY);

    Zone2 *pZoneGlobal(dt.createZone(NULL, ZL_GLOBAL));

    Zone2 *pBoundedZone(pZoneGlobal->convertToBoundedZone());
    dt.refine(pBoundedZone, min_angle, min_edge_tri, max_edge_tri, false);

    YourPeelPredicate_s decider;
    Zone2 *pPeeledZone = peelOffIf(pBoundedZone, &decider, false);
    if (pPeeledZone == NULL)
    {
        std::cout << "It seems a faulty predicate has peeled off all triangles, stop" << endl;
        exit(0);
    }

    std::vector<Triangle2 *> vTriangles;
    pPeeledZone->getTriangles(vTriangles);
    //std::cout << "after peelOff, size of vTriangles: " << vTriangles.size() << "\n";

    //std::vector<Point2 *> vVertices;
    //pPeeledZone->getVertices(vVertices);
    //std::cout << "after peelOff, size of vVertices: " << vVertices.size() << "\n";

    std::map<std::pair<double, double>, int> Map_Tri_JXY;
    //std::vector<RowVector6i> JM;
    //std::vector<Eigen::Vector3d> JM_linear;
    JM.resize(vTriangles.size());
    //JM_linear.resize(vTriangles.size());
    int pnt_ID = 0;

    for (size_t i = 0; i < vTriangles.size(); ++i)
    {
        RowVector6i tmp1;
        tmp1 << -1, -1, -1, -1, -1, -1;
        JM[i] = tmp1;
        for (size_t j = 0; j < 3; ++j)
        {
            Point2 *A = vTriangles[i]->getCorner(j); // the point coordinates
            std::pair<double, double> B;
            B.first = round(A->x(), 4);
            B.second = round(A->y(), 4);
            std::pair<std::pair<double, double>, int> SF = std::make_pair(B, pnt_ID);

            std::pair<std::map<std::pair<double, double>, int>::iterator, bool> ret = Map_Tri_JXY.insert(SF);

            size_t jk = 0;
            if (j == 1)
                jk = 2;
            else if (j == 2)
                jk = 4;

            if (ret.second)
            {
                //JM_linear[i](j) = pnt_ID;
                JM[i](0, jk) = pnt_ID;
                // std::cout << "add PNT ID = " << pnt_ID << std::endl;
                pnt_ID++;
            }
            else
            {
                std::map<std::pair<double, double>, int>::iterator it_s;
                it_s = Map_Tri_JXY.find(B);
                JM[i](0, jk) = it_s->second;
                //JM_linear[i](j) = it_s->second;
            }
        }
    }

    /*
    std::cout << "JM:\n";
    for (size_t i = 0; i < JM.size(); ++i)
    {

        for (size_t j = 0; j < 3; ++j)
        {
            size_t jk = 0;
            if (j == 1)
                jk = 2;
            else if (j == 2)
                jk = 4;
            std::cout << JM[i](0, jk) << ", ";
        }
        std::cout << "; \n";
    }
    */
    //std::vector<Eigen::Vector2d> JXY;
    JXY.resize(Map_Tri_JXY.size());
    std::map<std::pair<double, double>, int>::iterator it = Map_Tri_JXY.begin();
    while (it != Map_Tri_JXY.end())
    {
        Eigen::Vector2d UO;
        UO(0) = it->first.first;
        UO(1) = it->first.second;
        JXY[it->second] = UO;
        it++;
    }
    NO_Nodes_p = JXY.size();
    Find_neighbor_ele_and_shared_edge(JM, neigh_shared);
    Insert_central_pnt(JXY, JM, neigh_shared);
};

inline void Mesh_Polygon_2D_with_traces::Find_neighbor_ele_and_shared_edge(const std::vector<RowVector6i> JM, std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared)
{
    neigh_shared.resize(JM.size());
    for (size_t i = 0; i < JM.size(); ++i)
    {
        size_t A0, B0, C0;
        A0 = JM[i](0, 0);
        B0 = JM[i](0, 2);
        C0 = JM[i](0, 4);

        size_t A1, B1, C1;
        for (size_t j = 0; j < JM.size(); ++j)
        {
            if (j != i)
            {
                std::vector<size_t> VUY1, VUY0;
                A1 = JM[j](0, 0);
                B1 = JM[j](0, 2);
                C1 = JM[j](0, 4);
                /*
                std::cout << A0 << ", " << B0 << ", " << C0 << std::endl;
                std::cout << A1 << ", " << B1 << ", " << C1 << std::endl;
                std::cout << ";;;;;;\n";
                */
                size_t tag2 = 0;
                if (A0 == A1)
                {
                    tag2++;
                    VUY0.push_back(0);
                    VUY1.push_back(0);
                }
                else if (A0 == B1)
                {
                    tag2++;
                    VUY0.push_back(0);
                    VUY1.push_back(1);
                }
                else if (A0 == C1)
                {
                    tag2++;
                    VUY0.push_back(0);
                    VUY1.push_back(2);
                }

                if (B0 == A1)
                {
                    tag2++;
                    VUY0.push_back(1);
                    VUY1.push_back(0);
                }
                else if (B0 == B1)
                {
                    tag2++;
                    VUY0.push_back(1);
                    VUY1.push_back(1);
                }
                else if (B0 == C1)
                {
                    tag2++;
                    VUY0.push_back(1);
                    VUY1.push_back(2);
                }

                if (C0 == A1)
                {
                    tag2++;
                    VUY0.push_back(2);
                    VUY1.push_back(0);
                }
                else if (C0 == B1)
                {
                    tag2++;
                    VUY0.push_back(2);
                    VUY1.push_back(1);
                }
                else if (C0 == C1)
                {
                    tag2++;
                    VUY0.push_back(2);
                    VUY1.push_back(2);
                }
                /*
                if (tag2 == 2 && (VUY0.size() != 2 || VUY1.size() != 2))
                {
                    std::cout << "VUY0 size: " << VUY0.size() << std::endl;
                    std::cout << "VUY1 size: " << VUY1.size() << std::endl;
                    exit(0);
                }*/

                if (tag2 == 2 && VUY0.size() == 2 && VUY1.size() == 2)
                {
                    std::pair<size_t, Eigen::Vector2d> PUK;
                    PUK.first = j;
                    int Edge_no = -1;
                    if ((VUY0[0] == 0 && VUY0[1] == 1) || (VUY0[0] == 1 && VUY0[1] == 0))
                    {
                        Edge_no = 0;
                    }
                    else if ((VUY0[0] == 1 && VUY0[1] == 2) || (VUY0[0] == 2 && VUY0[1] == 1))
                    {
                        Edge_no = 1;
                    }
                    else if ((VUY0[0] == 0 && VUY0[1] == 2) || (VUY0[0] == 2 && VUY0[1] == 0))
                    {
                        Edge_no = 2;
                    }
                    PUK.second(0) = Edge_no;

                    //-----------
                    Edge_no = -1;
                    if ((VUY1[0] == 0 && VUY1[1] == 1) || (VUY1[0] == 1 && VUY1[1] == 0))
                    {
                        Edge_no = 0;
                    }
                    else if ((VUY1[0] == 1 && VUY1[1] == 2) || (VUY1[0] == 2 && VUY1[1] == 1))
                    {
                        Edge_no = 1;
                    }
                    else if ((VUY1[0] == 0 && VUY1[1] == 2) || (VUY1[0] == 2 && VUY1[1] == 0))
                    {
                        Edge_no = 2;
                    }
                    PUK.second(1) = Edge_no;

                    neigh_shared[i].push_back(PUK);
                }
                else if (tag2 > 2)
                {
                    std::cout << "Class: Mesh_Polygon_2D_with_traces. Error! shared edge just have two vertices!\n";
                    cout << "1st ele: " << A0 << ", " << B0 << ", "<< C0 << endl;
                    cout << "2nd ele: " << A1 << ", " << B1 << ", "<< C1 << endl;
                    cout << "tag2: " << tag2 << endl;
                    exit(0);
                }
            }
        }
    }

    // std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared
    /*
        for (size_t i = 0; i < neigh_shared.size(); ++i)
        {
            std::cout << "element NO " << i << ";\n";
            for (size_t j = 0; j < neigh_shared[i].size(); ++j)
            {
                std::cout << "\tele: " << neigh_shared[i][j].first << "; shared edged: " << neigh_shared[i][j].second(0) << ", " << neigh_shared[i][j].second(1) << std::endl;
            }
        }*/
};

inline void Mesh_Polygon_2D_with_traces::Insert_central_pnt(std::vector<Eigen::Vector2d> &JXY, std::vector<RowVector6i> &JM, const std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> neigh_shared)
{
    int update_pnt_ID = JXY.size();
    for (size_t i = 0; i < JM.size(); ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            size_t jk = 0;
            if (j == 0)
                jk = 1;
            else if (j == 1)
                jk = 3;
            else if (j == 2)
                jk = 5;

            if (JM[i](0, jk) == -1)
            {
                Eigen::Vector2d FG;
                FG(0) = .5 * (JXY[JM[i](0, jk - 1)](0) + JXY[JM[i](0, (jk + 1) % 6)](0));
                FG(1) = .5 * (JXY[JM[i](0, jk - 1)](1) + JXY[JM[i](0, (jk + 1) % 6)](1));
                JM[i](0, jk) = update_pnt_ID;
                size_t edge1;
                if (jk == 1)
                    edge1 = 0;
                else if (jk == 3)
                    edge1 = 1;
                else if (jk == 5)
                    edge1 = 2;
                //size_t ele1 = i;
                //std::cout << "now add pnt to elemen " << i << ", edge is " << edge1 << ", pnt No is: " << jk << ", define as " << JM[i](0, jk) << std::endl;
                for (size_t k = 0; k < neigh_shared[i].size(); ++k)
                {
                    if (edge1 == neigh_shared[i][k].second(0))
                    {
                        size_t ele2 = neigh_shared[i][k].first;
                        size_t edg2 = neigh_shared[i][k].second(1);
                        /*
                        std::cout << "\tele1: " << ele1 << "; "
                                  << "edge1: " << edge1 << "; ";
                        std::cout << "ele2: " << ele2 << "; "
                                  << "edge2: " << edg2 << "\n";*/
                        size_t sorted_node = 0;
                        if (edg2 == 0)
                            sorted_node = 1;
                        else if (edg2 == 1)
                            sorted_node = 3;
                        else if (edg2 == 2)
                            sorted_node = 5;

                        if (JM[ele2](0, sorted_node) == -1)
                            JM[ele2](0, sorted_node) = update_pnt_ID;
                        else
                        {
                            std::cout << "ERROR! the same pnt should have only one ID!\n";
                            exit(0);
                        }
                    }
                }
                JXY.push_back(FG);
                update_pnt_ID++;
            }
        }
    }
};

inline void Mesh_Polygon_2D_with_traces::Identifying_Frac_bound(const DFN::Polygon_convex_2D polygon, const vector<DFN::Line_seg_2D> Traces)
{
    Pnt_attribute.resize(JXY.size());

    for (size_t i = 0; i < JXY.size(); ++i)
    {
        DFN::Point_2D ThisPnt{JXY[i]};

        // polygon
        for (size_t k = 0; k < polygon.Corners.size(); ++k)
        {
            std::vector<Vector2d> Line_seg(2);
            Line_seg[0] << polygon.Corners[k](0), polygon.Corners[k](1);
            Line_seg[1] << polygon.Corners[(k + 1) % polygon.Corners.size()](0), polygon.Corners[(k + 1) % polygon.Corners.size()](1);
            if (ThisPnt.If_lies_on_a_line_seg(Line_seg) == true)
            {
                Pnt_attribute[i].If_frac_bound = true;
                break;
            };
        }

        // trace
        for (size_t j = 0; j < Traces.size(); ++j)
        {
            std::vector<Vector2d> Verts_2(2);
            Verts_2[0] << Traces[j].Point[0](0), Traces[j].Point[0](1);
            Verts_2[1] << Traces[j].Point[1](0), Traces[j].Point[1](1);
            if (ThisPnt.If_lies_on_a_line_seg(Verts_2) == true)
            {
                Pnt_attribute[i].If_trace = true;
                Trace_Tag.push_back(i);
                break;
            }
        }
    }
};

inline void Mesh_Polygon_2D_with_traces::Matlab_plot(string FileKey, const DFN::Splitting_Polygon_convex_2D_with_traces Poly_t)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "JXY=[";
    for (size_t j = 0; j < JXY.size(); ++j)
    {
        oss << JXY[j](0) << ", ";
        oss << JXY[j](1) << "; ";
    }
    oss << "];\n";

    oss << "JM";
    oss << "=[";
    for (size_t j = 0; j < JM.size(); ++j)
    {
        oss << JM[j](0) + 1 << ", ";
        oss << JM[j](1) + 1 << ", ";
        oss << JM[j](2) + 1 << ", ";
        oss << JM[j](3) + 1 << ", ";
        oss << JM[j](4) + 1 << ", ";
        oss << JM[j](5) + 1 << "; ";
    }
    oss << "];\n";

    oss << "Data = zeros(" << JXY.size() << ", 1);\n";

    oss << "figure(1)\n";
    oss << "P = patch('Vertices', JXY, 'Faces', JM, 'FaceVertexCData', Data, 'FaceColor', 'interp', 'EdgeAlpha', 0.2, 'facealpha', 0);\n\n\n";
    oss << "hold on;\n";
    oss << "Pnt_sets2 = [";
    for (size_t i = 0; i < Poly_t.Pnt_sets.size(); ++i)
    {
        oss << Poly_t.Pnt_sets[i](0) << ", " << Poly_t.Pnt_sets[i](1) << "; ";
    }
    oss << "];\n";
    oss << "hold on;\n";

    for (size_t i = 0; i < Poly_t.Topo_struct.size(); ++i)
    {
        oss << "plot([Pnt_sets2(" << Poly_t.Topo_struct[i].first + 1 << ",1), Pnt_sets2(" << Poly_t.Topo_struct[i].second + 1 << ",1)],[Pnt_sets2(" << Poly_t.Topo_struct[i].first + 1 << ",2), Pnt_sets2(" << Poly_t.Topo_struct[i].second + 1 << ",2)], 'LineWidth',2);\n";
        oss << "hold on;\n";
    }

    for (size_t i = 0; i < Pnt_attribute.size(); ++i)
    {
        if (Pnt_attribute[i].If_frac_bound == true)
        {

            oss << "text(" << JXY[i](0) << ", " << JXY[i][1] << ", 'x');\nhold on;\n";
        }

        if (Pnt_attribute[i].If_trace == true)
        {

            oss << "text(" << JXY[i](0) << ", " << JXY[i][1] << ", 't');\nhold on;\n";
        }
    }
    oss.close();
};

}; // namespace DFN