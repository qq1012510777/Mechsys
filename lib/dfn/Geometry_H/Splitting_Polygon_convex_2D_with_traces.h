#pragma once

#include "../Quaternion_H/Quaternion.h"
#include "Eigen/Dense"
#include "Polygon_convex_2D_with_traces.h"
#include <cmath>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Splitting_Polygon_convex_2D_with_traces
{
public:
    std::vector<Vector2d> Pnt_sets;
    std::vector<pair<size_t, size_t>> Topo_struct;
    std::vector<If_trace_segment> Topo_struct_attri;
    std::vector<Vector2d> Trace_pnts;

public:
    Splitting_Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D_with_traces Poly_trace, const double sub_edge, const std::vector<DFN::Line_seg_2D> Traces_Line);
    void Remove_isolated_and_tiny_constraints();
    void Matlab_plot(string FileKey);
};

inline Splitting_Polygon_convex_2D_with_traces::Splitting_Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D_with_traces Poly_trace, const double sub_edge, const std::vector<DFN::Line_seg_2D> Traces_Line)
{

    for (size_t i = 0; i < Traces_Line.size(); ++i)
    {
        if ((Traces_Line[i].Point[0] - Traces_Line[i].Point[1]).norm() < 0.01)
            Trace_pnts.push_back(Traces_Line[i].Point[0]);
    }

    this->Pnt_sets.resize(Poly_trace.Pnt_sets.size());
    for (size_t i = 0; i < this->Pnt_sets.size(); ++i)
        Pnt_sets[i] = Poly_trace.Pnt_sets[i].Coordinate;

    for (size_t i = 0; i < Poly_trace.Topo_struct.size(); ++i)
    {
        size_t pnt_A = Poly_trace.Topo_struct[i].first, pnt_B = Poly_trace.Topo_struct[i].second;

        if (pnt_A != pnt_B)
        {
            if (Poly_trace.Topo_struct_attri[i] == true)
            {
                DFN::Line_seg_2D this_line{Poly_trace.Pnt_sets[pnt_A].Coordinate, Poly_trace.Pnt_sets[pnt_B].Coordinate};

                vector<Vector2d> extra_pnts = this_line.Exact_pnts_along_this_line_seg_evenly(sub_edge);

                //std::vector<pair<size_t, size_t>> line_seg(1 + extra_pnts.size());

                if (1 + extra_pnts.size() == 1)
                {
                    Trace_pnts.push_back(Poly_trace.Pnt_sets[pnt_A].Coordinate);
                    Trace_pnts.push_back(Poly_trace.Pnt_sets[pnt_B].Coordinate);

                    Topo_struct.push_back(make_pair(pnt_A, pnt_B));
                    Topo_struct_attri.push_back(true);
                }
                else
                {
                    Pnt_sets.push_back(extra_pnts[0]);
                    Topo_struct.push_back(make_pair(pnt_A, Pnt_sets.size() - 1));
                    Topo_struct_attri.push_back(true);

                    if (extra_pnts.size() > 1)
                    {
                        //cout << "1" << endl;
                        for (size_t j = 0; j < extra_pnts.size() - 1; ++j)
                        {
                            //cout << "1" << endl;
                            Pnt_sets.push_back(extra_pnts[j + 1]);
                            Topo_struct.push_back(make_pair(Pnt_sets.size() - 2, Pnt_sets.size() - 1));
                            Topo_struct_attri.push_back(true);
                        }
                    }

                    //Pnt_sets.push_back(extra_pnts[extra_pnts.size() - 1]);
                    Topo_struct.push_back(make_pair(Pnt_sets.size() - 1, pnt_B));
                    Topo_struct_attri.push_back(true);

                    //
                    Trace_pnts.push_back(Poly_trace.Pnt_sets[pnt_A].Coordinate);
                    Trace_pnts.push_back(Poly_trace.Pnt_sets[pnt_B].Coordinate);
                    Trace_pnts.insert(Trace_pnts.end(), extra_pnts.begin(), extra_pnts.end());
                    // and their mid point
                    /*
                Trace_pnts.push_back((Poly_trace.Pnt_sets[pnt_A].Coordinate + extra_pnts[0]) / 2);
                Trace_pnts.push_back((Poly_trace.Pnt_sets[pnt_A].Coordinate + extra_pnts[extra_pnts.size() - 1]) / 2);

                for (size_t j = 1; j < extra_pnts.size(); ++j)
                {
                    Trace_pnts.push_back((extra_pnts[j] + extra_pnts[j - 1]) / 2);   
                }*/
                }
            }
            else
            {
                //cout << "polygon constraint\n";
                DFN::Line_seg_2D this_line{Poly_trace.Pnt_sets[pnt_A].Coordinate, Poly_trace.Pnt_sets[pnt_B].Coordinate};

                vector<Vector2d> extra_pnts = this_line.Exact_pnts_along_this_line_seg_evenly(1e10 /*sub_polygon*/);

                //std::vector<pair<size_t, size_t>> line_seg(1 + extra_pnts.size());

                if (1 + extra_pnts.size() == 1)
                {
                    Topo_struct.push_back(make_pair(pnt_A, pnt_B));
                    Topo_struct_attri.push_back(false);
                }
                else
                {
                    Pnt_sets.push_back(extra_pnts[0]);
                    Topo_struct.push_back(make_pair(pnt_A, Pnt_sets.size() - 1));
                    Topo_struct_attri.push_back(false);

                    if (extra_pnts.size() > 1)
                    {
                        //cout << "1" << endl;
                        for (size_t j = 0; j < extra_pnts.size() - 1; ++j)
                        {
                            //cout << "1" << endl;
                            Pnt_sets.push_back(extra_pnts[j + 1]);
                            Topo_struct.push_back(make_pair(Pnt_sets.size() - 2, Pnt_sets.size() - 1));
                            Topo_struct_attri.push_back(false);
                        }
                    }

                    //Pnt_sets.push_back(extra_pnts[extra_pnts.size() - 1]);
                    Topo_struct.push_back(make_pair(Pnt_sets.size() - 1, pnt_B));
                    Topo_struct_attri.push_back(false);
                }
            }
        }
    }

    //-------------------------------------------------------------------------
    for (size_t i = 0; i < Poly_trace.Topo_struct.size(); ++i)
    {
        size_t pnt_A = Poly_trace.Topo_struct[i].first, pnt_B = Poly_trace.Topo_struct[i].second;

        if (pnt_A == pnt_B && Poly_trace.Topo_struct_attri[i] == true)
        {
            Pnt_sets.push_back(Poly_trace.Pnt_sets[pnt_A].Coordinate);
            Topo_struct.push_back(make_pair(Pnt_sets.size() - 1, Pnt_sets.size() - 1));
            Topo_struct_attri.push_back(true);
        }
    }
    //------------------------------------------------------------------------
};

void Splitting_Polygon_convex_2D_with_traces::Matlab_plot(string FileKey)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";

    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        size_t pntID1 = Topo_struct[i].first, pntID2 = Topo_struct[i].second;

        if (Topo_struct_attri[i] == false)
        {
            oss << "p" << i << " = plot([" << round(Pnt_sets[pntID1](0), 4) << ", " << round(Pnt_sets[pntID2](0), 4) << "], [" << round(Pnt_sets[pntID1](1), 4) << ", " << round(Pnt_sets[pntID2](1), 4) << "], 'k-x', 'linewidth', 4);\n";
            oss << "hold on;\n";
        }
        else
        {
            oss << "p" << i << " = plot([" << round(Pnt_sets[pntID1](0), 4) << ", " << round(Pnt_sets[pntID2](0), 4) << "], [" << round(Pnt_sets[pntID1](1), 4) << ", " << round(Pnt_sets[pntID2](1), 4) << "], 'linewidth', 2);\n";
            oss << "hold on;\n";
        }
    }

    oss.close();
};

}; // namespace DFN