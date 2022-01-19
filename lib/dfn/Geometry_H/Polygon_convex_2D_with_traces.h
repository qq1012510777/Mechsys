#pragma once

#include "Eigen/Dense"
#include "Distance_2D.h"
#include "Line_seg_2D.h"
#include "Point_2D.h"
#include "Polygon_convex_2D.h"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

typedef bool If_trace_segment;

bool compVector2d(Vector2d &a, Vector2d &b)
{
    return a.norm() < b.norm();
}

namespace DFN
{

class Polygon_convex_2D_with_traces
{
public:
    std::vector<DFN::Point_2D> Pnt_sets;           // each Point_2D has a NO. which is the IDX of this vector
    std::vector<pair<size_t, size_t>> Topo_struct; // each pair is a segment from pnt NO. n1 to n2; if segment is a pnt actually, then n1 = n2
    std::vector<If_trace_segment> Topo_struct_attri;
    size_t Tag_for_Polygon_Pnt;

public:
    Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D polygon,
                                  const vector<DFN::Line_seg_2D> Traces);
    void Polygon_structure(const DFN::Polygon_convex_2D polygon);

    //splitting the polygon
    void Insert_IntersectionPnt_between_edge_and_trace(const vector<DFN::Line_seg_2D> Traces);

    void Split_traces(const vector<DFN::Line_seg_2D> Traces, vector<DFN::Line_seg_2D> &Split_Traces);

    size_t Insert_a_pnt_to_Pnt_sets(Vector2d A);

    void Insert_a_splitting_trace_to_Topo_struct(const size_t ID1, const size_t ID2);

    void Insert_a_splitting_trace(const DFN::Polygon_convex_2D polygon, DFN::Line_seg_2D splitting_trace);

    void Set_Topo_struct_attri(const DFN::Polygon_convex_2D polygon, const vector<DFN::Line_seg_2D> Traces);

    void Pointliked_traces(const vector<DFN::Line_seg_2D> Traces, const DFN::Polygon_convex_2D polygon); // point_like traces and do not intersect any traces or polygon edges

    void Matlab_plot(string FileKey);
};

inline Polygon_convex_2D_with_traces::Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D polygon,
                                                                    const vector<DFN::Line_seg_2D> Traces)
{
    //// create polygon structure first
    this->Polygon_structure(polygon);

    this->Insert_IntersectionPnt_between_edge_and_trace(Traces);
    Tag_for_Polygon_Pnt = Pnt_sets.size() - 1;
    //
    std::vector<DFN::Line_seg_2D> Split_Traces;
    this->Split_traces(Traces, Split_Traces);

    // insert trace to structure
    for (size_t i = 0; i < Split_Traces.size(); ++i)
    {
        this->Insert_a_splitting_trace(polygon, Split_Traces[i]);
    };

    this->Set_Topo_struct_attri(polygon, Traces);
    this->Pointliked_traces(Traces, polygon);
};

inline void Polygon_convex_2D_with_traces::Polygon_structure(const DFN::Polygon_convex_2D polygon)
{
    Pnt_sets.resize(polygon.Corners.size());

    for (size_t i = 0; i < polygon.Corners.size(); ++i)
    {
        Pnt_sets[i].Re_constructor(polygon.Corners[i]);
    }

    // remove overlapping point
    for (size_t i = 0; i < polygon.Corners.size();)
    {
        DFN::Line_seg_2D LineE{Pnt_sets[i].Coordinate,
                               Pnt_sets[(i + 1) % polygon.Corners.size()].Coordinate};

        if (LineE.If_is_a_point() == true)
        {
            Pnt_sets.erase(Pnt_sets.begin() + i + 1); // remove the 2nd point
        }
        else
            i++;
    }

    // structure
    Topo_struct.resize(Pnt_sets.size());
    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        Topo_struct[i] = std::make_pair(i, (i + 1) % Pnt_sets.size());
    }
};

inline void Polygon_convex_2D_with_traces::Insert_IntersectionPnt_between_edge_and_trace(const vector<DFN::Line_seg_2D> Traces)
{
    std::vector<Vector2d> Intersection_point;

    for (size_t i = 0; i < Traces.size(); ++i)
    {
        //DFN::Line_seg_2D Thisline_A{Traces[i].Point[0], Traces[i].Point[1]};

        DFN::Point_2D A1{Traces[i].Point[0]};
        DFN::Point_2D A2{Traces[i].Point[1]};

        for (size_t j = 0; j < Topo_struct.size(); ++j)
        {
            std::vector<Vector2d> Intersection_A;
            DFN::Line_seg_2D Thisline_B{Pnt_sets[Topo_struct[j].first].Coordinate, Pnt_sets[Topo_struct[j].second].Coordinate};
            //bool pu = Thisline_A.Intersection_between_two_lines(Thisline_B, Intersection_A);

            bool uit1 = A1.If_lies_on_a_line_seg(Thisline_B.Point);
            bool uit2 = A2.If_lies_on_a_line_seg(Thisline_B.Point);

            if (uit1 == true)
                Intersection_A.push_back(Traces[i].Point[0]);
            if (uit2 == true)
                Intersection_A.push_back(Traces[i].Point[1]);

            //if (pu == true)
            for (size_t k = 0; k < Intersection_A.size(); ++k)
            {
                Intersection_point.push_back(Intersection_A[k]);
            }
        }
    }

    //------------------------
    //cout << Intersection_point.size() << endl;
    for (size_t tk = 0; tk < Intersection_point.size(); ++tk)
    {

        Vector2d Pnt = Intersection_point[tk];
        //cout << Intersection_point[tk].transpose() << endl;
        // if the point is on the polygon bound, then change the polygon structure

        // first, make sure if the point is repetitive?
        bool repetitive = false;
        for (size_t i = 0; i < Pnt_sets.size(); ++i)
        {
            Vector2d AS = Pnt - Pnt_sets[i].Coordinate;
            if (abs(AS.norm()) < 0.01)
            {
                repetitive = true;
                break;
            }
        }
        if (repetitive == true)
            continue;

        // then, if the point splits the polygon and trace, change the polygon structure
        DFN::Point_2D AS{Pnt};

        for (size_t i = 0; i < Topo_struct.size(); ++i)
        {
            std::vector<Vector2d> LineE = {Pnt_sets[Topo_struct[i].first].Coordinate,
                                           Pnt_sets[Topo_struct[i].second].Coordinate};

            if (AS.If_lies_on_a_line_seg(LineE) == true)
            {
                //cout << 11 << endl;
                Pnt_sets.push_back(AS); // the point might touch other traces

                size_t PntID = Pnt_sets.size() - 1;
                std::pair<size_t, size_t> OPS = std::make_pair(PntID, Topo_struct[i].second);

                Topo_struct[i].second = PntID;

                Topo_struct.insert(Topo_struct.begin() + i + 1, OPS);
                break;
            }
            else if (AS.If_lies_on_a_line_seg(LineE) == false && i == Topo_struct.size() - 1)
            {
                cout << "error in Class 'Polygon_convex_2D_with_traces', function 'Insert_IntersectionPnt_between_edge_and_trace'!\n";
                cout << "should lie on one edge!\n";
                throw Error_throw_ignore("error in Class 'Polygon_convex_2D_with_traces', function 'Insert_IntersectionPnt_between_edge_and_trace'!\n");
            }
        }
    }
};

inline void Polygon_convex_2D_with_traces::Split_traces(const vector<DFN::Line_seg_2D> Traces, vector<DFN::Line_seg_2D> &Split_Traces)
{
    for (size_t i = 0; i < Traces.size(); ++i)
    {
        std::vector<Vector2d> Split_seg(1);
        Split_seg[0] = Traces[i].Point[0];

        DFN::Line_seg_2D LineA{Traces[i].Point[0], Traces[i].Point[1]};
        //cout << "Trace " << i << ": " << Traces[i].Point[0].transpose() << ", " << Traces[i].Point[1].transpose() << endl;

        for (size_t j = 0 && j != i; j < Traces.size(); ++j)
        {
            DFN::Line_seg_2D LineB{Traces[j].Point[0], Traces[j].Point[1]};
            //cout << "\tTrace " << j << ": " << Traces[j].Point[0].transpose() << ", " << Traces[j].Point[1].transpose() << endl;
            std::vector<Vector2d> Intersection;
            bool pot = LineA.Intersection_between_two_lines(LineB, Intersection);

            if (pot == true)
                for (size_t k = 0; k < Intersection.size(); ++k)
                {
                    Split_seg.push_back(Intersection[k]);
                    //cout << "\t\t" << Intersection[k].transpose() << endl;
                }
        }
        Split_seg.push_back(Traces[i].Point[1]);

        //sort
        Vector2d XU = Split_seg[0];
        for (size_t j = 0; j < Split_seg.size(); ++j)
            Split_seg[j] = Split_seg[j] - XU;

        sort(Split_seg.begin(), Split_seg.end(), compVector2d);

        for (size_t j = 0; j < Split_seg.size(); ++j)
            Split_seg[j] = Split_seg[j] + XU;
        //sort finished

        //remove overlapping point
        for (size_t j = 0; j < Split_seg.size() - 1;)
        {
            Vector2d AY = Split_seg[j] - Split_seg[j + 1];

            if (AY.norm() < 0.01)
            {
                Split_seg.erase(Split_seg.begin() + j + 1);
            }
            else
                j++;
        }

        // new
        for (size_t j = 0; j < Split_seg.size() - 1; ++j)
        {
            DFN::Line_seg_2D LineO{Split_seg[j], Split_seg[j + 1]};
            Split_Traces.push_back(LineO);
        }
    }
};

inline size_t Polygon_convex_2D_with_traces::Insert_a_pnt_to_Pnt_sets(Vector2d A)
{
    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        Vector2d YU = A - Pnt_sets[i].Coordinate;
        if (YU.norm() < 0.01)
        {
            return i;
        }
    }

    DFN::Point_2D PO{A};
    Pnt_sets.push_back(PO);
    return (Pnt_sets.size() - 1);
};

inline void Polygon_convex_2D_with_traces::Insert_a_splitting_trace_to_Topo_struct(const size_t ID1, const size_t ID2)
{
    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        if ((Topo_struct[i].first == ID1 && Topo_struct[i].second == ID2) ||
            (Topo_struct[i].first == ID2 && Topo_struct[i].second == ID1))
        {
            //existed
            return;
        }
    }
    Topo_struct.push_back(std::make_pair(ID1, ID2));
};

inline void Polygon_convex_2D_with_traces::Insert_a_splitting_trace(const DFN::Polygon_convex_2D polygon, DFN::Line_seg_2D splitting_trace)
{
    DFN::Line_seg_2D Trace_C{splitting_trace.Point[0], splitting_trace.Point[1]};

    for (size_t i = 0; i < polygon.Corners.size(); ++i)
    {
        std::vector<Vector2d> INterO;
        DFN::Line_seg_2D Trace_D{polygon.Corners[i], polygon.Corners[(i + 1) % polygon.Corners.size()]};
        bool POF = Trace_C.If_two_lines_overlap(Trace_D, INterO);
        if (POF == true)
        {
            return; // because the Trace_C overlaps with one edge of the polyon
        }
    }

    // if not overlap
    size_t ID1 = this->Insert_a_pnt_to_Pnt_sets(splitting_trace.Point[0]);
    size_t ID2 = this->Insert_a_pnt_to_Pnt_sets(splitting_trace.Point[1]);
    this->Insert_a_splitting_trace_to_Topo_struct(ID1, ID2);
};

inline void Polygon_convex_2D_with_traces::Set_Topo_struct_attri(const DFN::Polygon_convex_2D polygon, const vector<DFN::Line_seg_2D> Traces)
{
    Topo_struct_attri.resize(Topo_struct.size());

    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        Topo_struct_attri[i] = true;
        size_t pnt1 = Topo_struct[i].first, pnt2 = Topo_struct[i].second;
        if (pnt1 > Tag_for_Polygon_Pnt || pnt2 > Tag_for_Polygon_Pnt)
        {
            Topo_struct_attri[i] = true;
            continue;
        }

        DFN::Line_seg_2D Thisline_A{Pnt_sets[pnt1].Coordinate, Pnt_sets[pnt2].Coordinate};
        DFN::Point_2D MidPnt{(Pnt_sets[pnt1].Coordinate + Pnt_sets[pnt2].Coordinate) * 0.50};

        for (size_t j = 0; j < polygon.Corners.size(); ++j)
        {
            std::vector<Vector2d> Thisline = {polygon.Corners[j], polygon.Corners[(j + 1) % polygon.Corners.size()]};

            bool iflieon = MidPnt.If_lies_on_a_line_seg(Thisline);

            if (iflieon == true)
            {
                Topo_struct_attri[i] = false;
                break;
            }
        }

        for (size_t j = 0; j < Traces.size(); ++j)
        {
            std::vector<Vector2d> Thisline = Traces[j].Point;

            bool iflieon = MidPnt.If_lies_on_a_line_seg(Thisline);

            if (iflieon == true)
            {
                Topo_struct_attri[i] = true;
                break;
            }
        }
    }
};

inline void Polygon_convex_2D_with_traces::Pointliked_traces(const vector<DFN::Line_seg_2D> Traces, const DFN::Polygon_convex_2D polygon)
{
    for (size_t i = 0; i < Traces.size(); ++i)
    {
        DFN::Line_seg_2D Traces_this{Traces[i].Point[0], Traces[i].Point[1]};

        if (Traces_this.If_is_a_point() == true)
        {
            DFN::Point_2D Thispnt{Traces[i].Point[0]};
            bool IfLiesOn = false;
            for (size_t j = 0; j < polygon.Corners.size(); ++j)
            {
                DFN::Line_seg_2D SEGA{polygon.Corners[j], polygon.Corners[(j + 1) % polygon.Corners.size()]};
                bool yui = Thispnt.If_lies_on_a_line_seg(SEGA.Point);

                if (yui == true && SEGA.If_is_a_point() == false)
                {
                    IfLiesOn = true;
                    break;
                }
            }

            if (IfLiesOn == true)
            {
                continue;
            }

            for (size_t j = 0; j < Traces.size(); ++j)
            {
                if (j != i)
                {
                    DFN::Line_seg_2D SEGA{Traces[j].Point[0], Traces[j].Point[1]};
                    bool yui = Thispnt.If_lies_on_a_line_seg(SEGA.Point);

                    if (yui == true && SEGA.If_is_a_point() == false)
                    {
                        IfLiesOn = true;
                        break;
                    }
                }
            }

            if (IfLiesOn == false)
            {
                Pnt_sets.push_back(Thispnt);  
                Topo_struct.push_back(std::make_pair(Pnt_sets.size() - 1, Pnt_sets.size() - 1));
                Topo_struct_attri.push_back(true);
            }
        }
    }
};

void Polygon_convex_2D_with_traces::Matlab_plot(string FileKey)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "Pnt_sets = [";
    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        oss << Pnt_sets[i].Coordinate(0) << ", " << Pnt_sets[i].Coordinate(1) << "; ";
    }
    oss << "];\n";
    oss << "hold on;\n";

    //std::vector<string> IRE = {"-o", "-+", "-x"};
    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        if (Topo_struct_attri[i] == false)
        {
            oss << "P" << i + 1 << " = plot([Pnt_sets(" << Topo_struct[i].first + 1 << ",1), Pnt_sets(" << Topo_struct[i].second + 1 << ",1)],[Pnt_sets(" << Topo_struct[i].first + 1 << ",2), Pnt_sets(" << Topo_struct[i].second + 1 << ",2)], 'k-o', 'LineWidth',2);\n";
            oss << "hold on;\n";
        }
        else
        {
            oss << "P" << i + 1 << " = plot([Pnt_sets(" << Topo_struct[i].first + 1 << ",1), Pnt_sets(" << Topo_struct[i].second + 1 << ",1)],[Pnt_sets(" << Topo_struct[i].first + 1 << ",2), Pnt_sets(" << Topo_struct[i].second + 1 << ",2)], 'LineWidth',2);\n";
            oss << "hold on;\n";
        }
    }

    oss.close();
};

}; // namespace DFN