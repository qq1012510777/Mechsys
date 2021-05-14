#pragma once

#include "../Quaternion_H/Quaternion.h"
#include "Dense"
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
  

public:
    Splitting_Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D_with_traces Poly_trace, const double sub_edge);
};

inline Splitting_Polygon_convex_2D_with_traces::Splitting_Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D_with_traces Poly_trace, const double sub_edge)
{
    this->Pnt_sets.resize(Poly_trace.Pnt_sets.size());
    for (size_t i = 0; i < this->Pnt_sets.size(); ++i)
        Pnt_sets[i] = Poly_trace.Pnt_sets[i].Coordinate;

    for (size_t i = 0; i < Poly_trace.Topo_struct.size(); ++i)
    {
        size_t pnt_A = Poly_trace.Topo_struct[i].first, pnt_B = Poly_trace.Topo_struct[i].second;

        DFN::Line_seg_2D this_line{Poly_trace.Pnt_sets[pnt_A].Coordinate, Poly_trace.Pnt_sets[pnt_B].Coordinate};
        vector<Vector2d> extra_pnts = this_line.Exact_pnts_along_this_line_seg_evenly(sub_edge);

        //std::vector<pair<size_t, size_t>> line_seg(1 + extra_pnts.size());

        if (1 + extra_pnts.size() == 1)
        {
            Topo_struct.push_back(make_pair(pnt_A, pnt_B));
        }
        else
        {
            Pnt_sets.push_back(extra_pnts[0]);
            Topo_struct.push_back(make_pair(pnt_A, Pnt_sets.size() - 1));

            if (extra_pnts.size() > 1)
            {
                //cout << "1" << endl;
                for (size_t j = 0; j < extra_pnts.size() - 1; ++j)
                {
                    //cout << "1" << endl;
                    Pnt_sets.push_back(extra_pnts[j + 1]);
                    Topo_struct.push_back(make_pair(Pnt_sets.size() - 2, Pnt_sets.size() - 1));
                }
            }

            //Pnt_sets.push_back(extra_pnts[extra_pnts.size() - 1]);
            Topo_struct.push_back(make_pair(Pnt_sets.size() - 1, pnt_B));
        }
    }
};

}; // namespace DFN