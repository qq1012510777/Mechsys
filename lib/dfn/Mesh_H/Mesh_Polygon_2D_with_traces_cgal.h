#pragma once

//CGAL
#define CGAL_MESH_2_OPTIMIZER_VERBOSE
//#define CGAL_MESH_2_OPTIMIZERS_DEBUG
//#define CGAL_MESH_2_SIZING_FIELD_USE_BARYCENTRIC_COORDINATES
#include "CGAL/Constrained_Delaunay_triangulation_2.h"
#include "CGAL/Delaunay_mesh_face_base_2.h"
#include "CGAL/Delaunay_mesh_size_criteria_2.h"
#include "CGAL/Delaunay_mesh_vertex_base_2.h"
#include "CGAL/Delaunay_mesher_2.h"
#include "CGAL/Delaunay_mesher_no_edge_refinement_2.h"
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/lloyd_optimize_mesh_2.h"
#include <iostream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
//typedef CGAL::Delaunay_mesher_no_edge_refinement_2<CDT, Criteria> Mesher;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point_CGAL;

//----------

#include "../Geometry_H/Splitting_Polygon_convex_2D_with_traces.h"
#include "../Math_WL_H/Math_WL.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

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
                                const std::vector<If_trace_segment> Topo_struct_attri,
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
                                                                const std::vector<If_trace_segment> Topo_struct_attri,
                                                                const double min_angle,
                                                                const double min_edge_tri,
                                                                const double max_edge_tri,
                                                                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared,
                                                                size_t &NO_Nodes_p)
{
    std::vector<Point_CGAL> vInputPoints;
    vInputPoints.resize(Pnt_sets.size());

    //cout << "std::vector<Point> vInputPoints(" << Pnt_sets.size() << ");\n";
    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        vInputPoints[i] = Point_CGAL(Pnt_sets[i](0), Pnt_sets[i](1));
        //cout << "vInputPoints[" << i << "] = Point(" << Pnt_sets[i](0) << ", " << Pnt_sets[i](1) << ");\n";
    }
    //cout << "\n\n\n";
    CDT cdt;

    std::vector<Vertex_handle> va(128);
    for (size_t i = 0; i < va.size(); ++i)
        va[i] = cdt.insert(vInputPoints[i]);

    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        cdt.insert_constraint(va[Topo_struct[i].first], va[Topo_struct[i].second]);
        //cout << "cdt.insert_constraint(va[" <<Topo_struct[i].first << "], va[" << Topo_struct[i].second<<"]);\n";
    }
    cout << "meahing ...\n";
    Mesher mesher(cdt);
    cout << "1st meshing finished\n";
    cout << max_edge_tri << endl;
    mesher.set_criteria(Criteria(0.125, max_edge_tri));
    mesher.refine_mesh();
    cout << "2nd meshing finished\n";
    CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = 10);

    //-----------------

    std::map<std::pair<double, double>, int> Map_Tri_JXY;
    //std::vector<RowVector6i> JM;
    //std::vector<Eigen::Vector3d> JM_linear;
    JM.resize(cdt.number_of_faces());
    //JM_linear.resize(vTriangles.size());
    int pnt_ID = 0;

    CDT::Finite_faces_iterator f_iter;
    size_t iiio;
    for (f_iter = cdt.finite_faces_begin(), iiio = 0;
         f_iter != cdt.finite_faces_end();
         f_iter++, iiio++)
    {
        RowVector6i tmp1;
        tmp1 << -1, -1, -1, -1, -1, -1;
        JM[iiio] = tmp1;

        std::vector<size_t> PointID(3);

        for (size_t j = 0; j < 3; ++j)
        {
            Point_CGAL A = f_iter->vertex(j)->point(); // the point coordinates
            std::pair<double, double> B;
            B.first = round(A.x(), 4);
            B.second = round(A.y(), 4);

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
                JM[iiio](0, jk) = pnt_ID;
                PointID[j] = pnt_ID;
                // std::cout << "add PNT ID = " << pnt_ID << std::endl;
                pnt_ID++;
            }
            else
            {
                std::map<std::pair<double, double>, int>::iterator it_s;
                it_s = Map_Tri_JXY.find(B);
                JM[iiio](0, jk) = it_s->second;
                PointID[j] = it_s->second;
                //JM_linear[i](j) = it_s->second;
            }
        }

        if (PointID[0] == PointID[1] || PointID[0] == PointID[2] || PointID[1] == PointID[2])
        {
            cout << "Error! Class: Mesh_Polygon_2D_with_traces! find bad triangle!\n";
            Point_CGAL A0 = f_iter->vertex(0)->point();
            Point_CGAL A1 = f_iter->vertex(1)->point();
            Point_CGAL A2 = f_iter->vertex(2)->point();
            cout << "the first vertex: " << A0.x() << ", " << A0.y() << endl;
            cout << "the second vertex: " << A1.x() << ", " << A1.y() << endl;
            cout << "the third vertex: " << A2.x() << ", " << A2.y() << endl;

            cout << "The topo structure: ";
            cout << PointID[0] << ", " << PointID[1] << ", " << PointID[2] << endl;
            exit(0);
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
                    cout << "1st ele: " << A0 << ", " << B0 << ", " << C0 << endl;
                    cout << "2nd ele: " << A1 << ", " << B1 << ", " << C1 << endl;
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