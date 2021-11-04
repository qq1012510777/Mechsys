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

//
#include <gmsh.h>

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
    std::vector<size_t> Trace_Tag; // record trace point

    std::vector<bool> IfCorner;

public:
    Mesh_Polygon_2D_with_traces(const std::vector<Vector2d> Pnt_sets,
                                const std::vector<pair<size_t, size_t>> Topo_struct,
                                const std::vector<If_trace_segment> Topo_struct_attri,
                                const double subpolygon,
                                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared,
                                size_t &NO_Nodes_p,
                                DFN::Polygon_convex_2D polyframe,
                                DFN::Splitting_Polygon_convex_2D_with_traces splitted_polygon,
                                std::vector<DFN::Line_seg_2D> Traces_Line);

    void Matlab_plot(string FileKey, const DFN::Splitting_Polygon_convex_2D_with_traces Poly_t);

    void Identifying_Frac_bound(const DFN::Polygon_convex_2D polygon, const vector<DFN::Line_seg_2D> Traces);

    bool If_new_pnt_generates_and_splits_traces(std::vector<DFN::Line_seg_2D> Traces_Line, std::vector<Vector2d> Trace_pnts);

    void Identify_corner_middle_point(size_t &NO_Nodes_p);
};

inline Mesh_Polygon_2D_with_traces::Mesh_Polygon_2D_with_traces(const std::vector<Vector2d> Pnt_sets,
                                                                const std::vector<pair<size_t, size_t>> Topo_struct,
                                                                const std::vector<If_trace_segment> Topo_struct_attri,
                                                                const double subpolygon,
                                                                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> &neigh_shared,
                                                                size_t &NO_Nodes_p,
                                                                DFN::Polygon_convex_2D polyframe,
                                                                DFN::Splitting_Polygon_convex_2D_with_traces splitted_polygon,
                                                                std::vector<DFN::Line_seg_2D> Traces_Line)
{

    //----------------
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 2); // default level is 5
    gmsh::model::add("t2");

    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        std::vector<Vector2d>::iterator ity = find(splitted_polygon.Trace_pnts.begin(), splitted_polygon.Trace_pnts.end(), Vector2d{Pnt_sets[i](0), Pnt_sets[i](1)});

        if(ity != splitted_polygon.Trace_pnts.end())
            gmsh::model::geo::addPoint(Pnt_sets[i](0), Pnt_sets[i](1), 0, subpolygon, i + 1);
        else
            gmsh::model::geo::addPoint(Pnt_sets[i](0), Pnt_sets[i](1), 0, 0, i + 1);
    }

    size_t line_Tag = 1;

    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        gmsh::model::geo::addLine(Topo_struct[i].first + 1, Topo_struct[i].second + 1, i + 1);
        line_Tag = i + 1;
        if (Topo_struct[i].second == 0)
        {
            break;
        }
    }

    std::vector<int> curveloop(line_Tag);
    for (size_t i = 0; i < curveloop.size(); ++i)
        curveloop[i] = i + 1;
    gmsh::model::geo::addCurveLoop(curveloop, 1);

    std::vector<int> surfaceloop = {1};
    gmsh::model::geo::addPlaneSurface(surfaceloop, 1);

    gmsh::model::geo::synchronize();
    //cout << 1000 << endl;

    //------------------trace
    size_t SizeOfEmbededPoint = 0;
    for (size_t i = line_Tag; i < Topo_struct.size(); ++i)
    {
        if (Topo_struct[i].first != Topo_struct[i].second)
            gmsh::model::geo::addLine(Topo_struct[i].first + 1,
                                      Topo_struct[i].second + 1,
                                      i + 1);
        else
            SizeOfEmbededPoint++;
    }
    //cout << 1001 << endl;
    gmsh::model::geo::synchronize();

    std::vector<int> embeded_curveloop(Topo_struct.size() - line_Tag - SizeOfEmbededPoint);
    for (size_t i = line_Tag; i < Topo_struct.size() - SizeOfEmbededPoint; ++i)
    {
        embeded_curveloop[i - line_Tag] = i + 1;
        //cout << "embeded_curveloop[i] " << embeded_curveloop[i] << endl;
    }
    //cout << 1002 << endl;
    if (embeded_curveloop.size() > 0)
    {
        gmsh::model::mesh::embed(1, embeded_curveloop, 2, 1);
    }
    //cout << 1003 << endl;
    if (embeded_curveloop.size() > 0)
        gmsh::model::geo::synchronize();

    if (SizeOfEmbededPoint > 0)
    {
        gmsh::model::geo::synchronize();
        std::vector<int> embeded_pnt(SizeOfEmbededPoint);
        for (size_t i = 0; i < SizeOfEmbededPoint; ++i)
        {
            embeded_pnt[i] = Topo_struct[Topo_struct.size() - 1 - i].first + 1;
        }
        gmsh::model::mesh::embed(0, embeded_pnt, 2, 1);
        gmsh::model::geo::synchronize();
    }

    gmsh::model::mesh::setOrder(2);
    gmsh::option::setNumber("Mesh.ElementOrder", 2);

    gmsh::model::mesh::generate(2);

    std::vector<std::size_t> nodes;
    std::vector<double> coord, coordParam;
    gmsh::model::mesh::getNodes(nodes, coord, coordParam);

    this->JXY.resize(coord.size() / 3);
    for (size_t i = 0; i < this->JXY.size(); i++)
    {
        this->JXY[i] << coord[i * 3], coord[i * 3 + 1];
    }

    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 2, -1);

    this->JM.resize(elemNodeTags[0].size() / 6);

    for (size_t i = 0; i < this->JM.size(); i++)
    {
        JM[i] << elemNodeTags[0][i * 6 + 0] - 1,
            elemNodeTags[0][i * 6 + 3] - 1,
            elemNodeTags[0][i * 6 + 1] - 1,
            elemNodeTags[0][i * 6 + 4] - 1,
            elemNodeTags[0][i * 6 + 2] - 1,
            elemNodeTags[0][i * 6 + 5] - 1;
    }

    gmsh::clear();
    gmsh::finalize();

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

    //--------------
    this->Identifying_Frac_bound(polyframe, Traces_Line);
    this->Identify_corner_middle_point(NO_Nodes_p);

    bool If_new_pnt_splits_trace = this->If_new_pnt_generates_and_splits_traces(Traces_Line, splitted_polygon.Trace_pnts);

    if (If_new_pnt_splits_trace == true)
    {
        throw Error_throw_ignore("Error! new pnt generates and splits the trace!\n");  
    }
}

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

inline void Mesh_Polygon_2D_with_traces::Identify_corner_middle_point(size_t &NO_Nodes_p)
{
    NO_Nodes_p = 0;

    Eigen::VectorXd JXY_op;
    JXY_op = Eigen::VectorXd::Zero(this->JXY.size());
    IfCorner.resize(this->JXY.size());

    for (size_t i = 0; i < JM.size(); ++i)
    {
        for (size_t j = 0; j < 6; ++j)
        {
            size_t pntIU = JM[i][j];

            if (JXY_op[pntIU] == 0) // not accessed
            {
                JXY_op[pntIU] = 1;

                if (j % 2 != 0)
                {
                    IfCorner[pntIU] = false; // mid
                }
                else
                {
                    IfCorner[pntIU] = true; // corner
                    NO_Nodes_p++;
                };
            }
        }
    }
};

inline bool Mesh_Polygon_2D_with_traces::If_new_pnt_generates_and_splits_traces(std::vector<DFN::Line_seg_2D> Traces_Line, std::vector<Vector2d> Trace_pnts)
{
    for (size_t i = 0; i < this->JXY.size(); ++i)
    {
        if (this->IfCorner[i] == true && this->Pnt_attribute[i].If_trace == true)
        {
            DFN::Point_2D thisPNT{JXY[i]};

            for (size_t k = 0; k < Trace_pnts.size(); ++k)
            {
                Vector2d UIO = thisPNT.Coordinate - Trace_pnts[k];
                if ((abs(UIO(0)) < 0.05 && abs(UIO(1)) < 0.05) || UIO.norm() < 0.05)
                {
                    break; //
                }
                else
                {
                    if (k == Trace_pnts.size() - 1)
                    {
                        cout << "-----------------------\n";
                        cout << "\tNew point: ";
                        cout << thisPNT.Coordinate.transpose() << endl;

                        cout << "\nTrace point:\n";
                        for (size_t io = 0; io < Trace_pnts.size(); ++io)
                            cout << Trace_pnts[io].transpose() << endl;

                        cout << "\nTrace Line seg:\n";
                        for (size_t io = 0; io < Traces_Line.size(); ++io)
                            cout << Traces_Line[io].Point[0].transpose() << ", " << Traces_Line[io].Point[1].transpose() << endl;
                        cout << "-----------------------\n";

                        return true; // this pnt splits constraints
                    }
                }
            }
        }
    }
    return false;
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

    for (size_t i = 0; i < Poly_t.Trace_pnts.size(); ++i)
    {
        oss << "scatter(" << Poly_t.Trace_pnts[i](0) << ", " << Poly_t.Trace_pnts[i](1) << ", 'o');\n";
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