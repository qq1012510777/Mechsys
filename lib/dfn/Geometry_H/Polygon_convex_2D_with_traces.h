#pragma once

#include "Dense"
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

namespace DFN
{

class Polygon_convex_2D_with_traces
{
public:
    std::vector<Vector2d> Corners;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    std::vector<DFN::Point_2D> Pnt_sets;           // each Point_2D has a NO. which is the IDX of this vector
    std::vector<pair<size_t, size_t>> Topo_struct; // each pair is a segment from pnt NO. n1 to n2; if segment is a pnt actually, then n1 = n2

public:
    Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D polygon, const vector<DFN::Line_seg_2D> Traces);
    void Set_topo_struct(const vector<DFN::Line_seg_2D> Traces);
    bool Insert_a_pnt_into_topo_line(std::vector<pair<size_t, DFN::Point_2D>> &topoline, DFN::Point_2D pnt, size_t &idx);
    bool Check_repetitive_line_seg(const size_t A, const size_t B);
    bool Check_repetitive_pnt(const DFN::Point_2D G, size_t &idx);
    void Matlab_plot(string FileKey);
};

inline Polygon_convex_2D_with_traces::Polygon_convex_2D_with_traces(const DFN::Polygon_convex_2D polygon, const vector<DFN::Line_seg_2D> Traces)
{
    this->Corners = polygon.Corners;
    this->x_min = polygon.x_min;
    this->x_max = polygon.x_max;
    this->y_min = polygon.y_min;
    this->y_max = polygon.y_max;

    Pnt_sets.resize(polygon.Corners.size());

    for (size_t i = 0; i < this->Corners.size(); ++i)
    {
        Pnt_sets[i].Re_constructor(this->Corners[i]);
    }

    this->Set_topo_struct(Traces);
};

inline void Polygon_convex_2D_with_traces::Set_topo_struct(const vector<DFN::Line_seg_2D> Traces)
{
    std::vector<std::vector<pair<size_t, DFN::Point_2D>>> Polygon_topo_line_seg(Pnt_sets.size());

    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        Polygon_topo_line_seg[i].push_back(std::make_pair(i, DFN::Point_2D{Corners[i]}));
        Polygon_topo_line_seg[i].push_back(std::make_pair((i + 1) % Pnt_sets.size(), DFN::Point_2D{Corners[(i + 1) % Pnt_sets.size()]}));
    }

    std::vector<std::vector<pair<size_t, DFN::Point_2D>>> Trace_topo_line_seg(Traces.size());
    for (size_t i = 0; i < Traces.size(); ++i)
    {
        Trace_topo_line_seg[i].push_back(std::make_pair(1e10, DFN::Point_2D{Traces[i].Point[0]}));
        Trace_topo_line_seg[i].push_back(std::make_pair(1e10, DFN::Point_2D{Traces[i].Point[1]}));
    };
    for (size_t i = 0; i < this->Corners.size(); ++i)
    {
        DFN::Line_seg_2D edge_a{Corners[i], Corners[(i + 1) % Corners.size()]};

        for (size_t j = 0; j < Traces.size(); ++j)
        {
            DFN::Line_seg_2D trace_a{Traces[j].Point[0], Traces[j].Point[1]};

            if (trace_a.If_is_a_point() == false)
            {
                std::vector<Vector2d> Intersection_t_e;
                bool rt = edge_a.Intersection_between_two_lines(trace_a, Intersection_t_e);

                if (rt == true && Intersection_t_e.size() == 1)
                {
                    DFN::Point_2D this_pnt{Intersection_t_e[0]};
                    size_t iy = 0;
                    bool ut = this->Insert_a_pnt_into_topo_line(Polygon_topo_line_seg[i], this_pnt, iy);

                    if (ut == true)
                    {
                        //insert successfully
                        Pnt_sets.push_back(this_pnt);

                        Polygon_topo_line_seg[i][iy].first = Pnt_sets.size() - 1;
                        //how to handle trace
                        for (size_t k = 0; k < Trace_topo_line_seg[j].size(); ++k)
                        {
                            if (Trace_topo_line_seg[j][k].second.If_overlap_with_another_pnt(this_pnt) == true)
                            {
                                Trace_topo_line_seg[j][k].first = Pnt_sets.size() - 1;
                                break;
                            }
                        };
                    }
                    else if (ut == false && iy != 1e10)
                    {
                        //find intersection is overlapping with an end (iy)

                        for (size_t k = 0; k < Trace_topo_line_seg[j].size(); ++k)
                        {
                            if (Trace_topo_line_seg[j][k].second.If_overlap_with_another_pnt(this_pnt) == true)
                            {
                                Trace_topo_line_seg[j][k].first = iy;
                                break;
                            }
                        };
                    }
                    else if (ut == false && iy == 1e10)
                    {
                        cout << "Polygon_convex_2D_with_traces::Set_topo_struct, found an intersection between an edge of a polygon and a trace, but this intersection is not on the edge??\n";
                        cout << "type 1!\n";
                        exit(0);
                    }
                }
                else if (rt == true && Intersection_t_e.size() == 2)
                {
                    //intersection overlap with the edge
                    for (size_t k = 0; k < Intersection_t_e.size(); ++k)
                    {
                        DFN::Point_2D this_pnt{Intersection_t_e[k]};
                        size_t iy = 0;
                        bool ut = this->Insert_a_pnt_into_topo_line(Polygon_topo_line_seg[i], this_pnt, iy);

                        if (ut == true)
                        {
                            //insert successfully
                            Pnt_sets.push_back(this_pnt);

                            Polygon_topo_line_seg[i][iy].first = Pnt_sets.size() - 1;
                            //how to handle trace
                            for (size_t m = 0; m < Trace_topo_line_seg[j].size(); ++m)
                            {
                                if (Trace_topo_line_seg[j][m].second.If_overlap_with_another_pnt(this_pnt) == true)
                                {
                                    Trace_topo_line_seg[j][m].first = Pnt_sets.size() - 1;
                                    break;
                                }
                            };
                        }
                        else if (ut == false && iy != 1e10)
                        {
                            //find intersection is overlapping with an end

                            for (size_t m = 0; m < Trace_topo_line_seg[j].size(); ++m)
                            {
                                if (Trace_topo_line_seg[j][m].second.If_overlap_with_another_pnt(this_pnt) == true)
                                {
                                    Trace_topo_line_seg[j][m].first = iy;
                                    break;
                                }
                            };
                        }
                        else if (ut == false && iy == 1e10)
                        {
                            cout << "Polygon_convex_2D_with_traces::Set_topo_struct, found an intersection between an edge of a polygon and a trace, but this intersection is not on the edge??\n";
                            cout << "type 2!\n";
                            exit(0);
                        }
                    }
                }
                else if (rt == false)
                {
                    ;
                    // how to address trace
                }
            }
            else
            {
                Vector3d A, B;
                A << Corners[i](0), Corners[i](1), 0;
                B << Corners[(i + 1) % Corners.size()](0), Corners[(i + 1) % Corners.size()](1), 0;
                std::vector<Vector3d> Verts_1{A, B};
                DFN::Point_2D this_pnt{Traces[j].Point[0]};

                bool rw = this_pnt.If_lies_on_the_bounds_of_polygon(Verts_1);
                if (rw == true)
                {
                    size_t iy = 0;
                    bool ut = this->Insert_a_pnt_into_topo_line(Polygon_topo_line_seg[i], this_pnt, iy);

                    if (ut == true)
                    {
                        //insert successfully
                        Pnt_sets.push_back(this_pnt);

                        Polygon_topo_line_seg[i][iy].first = Pnt_sets.size() - 1;
                        //how to handle trace
                        for (size_t k = 0; k < Trace_topo_line_seg[j].size(); ++k)
                        {
                            if (Trace_topo_line_seg[j][k].second.If_overlap_with_another_pnt(this_pnt) == true)
                            {
                                Trace_topo_line_seg[j][k].first = Pnt_sets.size() - 1;
                                break;
                            }
                        };
                    }
                    else if (ut == false && iy != 1e10)
                    {
                        //find intersection is overlapping with an end
                        //how to handle trace

                        for (size_t k = 0; k < Trace_topo_line_seg[j].size(); ++k)
                        {
                            if (Trace_topo_line_seg[j][k].second.If_overlap_with_another_pnt(this_pnt) == true)
                            {
                                Trace_topo_line_seg[j][k].first = iy;
                                break;
                            }
                        };
                    }
                    else if (ut == false && iy == 1e10)
                    {
                        cout << "Polygon_convex_2D_with_traces::Set_topo_struct, found an intersection between an edge of a polygon and a trace, but this intersection is not on the edge??\n";
                        cout << "type 2!\n";
                        exit(0);
                    }
                }
                else
                {
                    //how to handle trace
                    ;
                }
            }
        }
    };

    for (size_t i = 0; i < Trace_topo_line_seg.size(); ++i)
    {

        for (size_t j = 0; j < Trace_topo_line_seg[i].size(); ++j)
        {
            if (Trace_topo_line_seg[i][j].first == 1e10)
            {
                Pnt_sets.push_back(Trace_topo_line_seg[i][j].second);

                Trace_topo_line_seg[i][j].first = Pnt_sets.size() - 1;
            }
        }
    }

    for (size_t i = 0; i < Polygon_topo_line_seg.size(); ++i)
    {
        for (size_t j = 0; j < Polygon_topo_line_seg[i].size() - 1; ++j)
        {
            Topo_struct.push_back(make_pair(Polygon_topo_line_seg[i][j].first, Polygon_topo_line_seg[i][j + 1].first));
        }
    }

    if (Trace_topo_line_seg.size() > 0)
        for (size_t i = 0; i < Trace_topo_line_seg.size() - 1; ++i)
        {

            DFN::Line_seg_2D A{Trace_topo_line_seg[i][0].second.Coordinate, Trace_topo_line_seg[i][Trace_topo_line_seg[i].size() - 1].second.Coordinate};

            //check if trace overlap with edge of polygon
            if (this->Check_repetitive_line_seg(Trace_topo_line_seg[i][0].first, Trace_topo_line_seg[i][1].first) == false)
            {
                // if this trace intersect other traces?
                for (size_t j = i + 1; j < Trace_topo_line_seg.size(); ++j)
                {
                    DFN::Line_seg_2D B{Trace_topo_line_seg[j][0].second.Coordinate, Trace_topo_line_seg[j][Trace_topo_line_seg[j].size() - 1].second.Coordinate};
                    std::vector<Vector2d> Intersection_AB;
                    bool cg = A.Intersection_between_two_lines(B, Intersection_AB);
                    if (cg == true && Intersection_AB.size() == 1)
                    {
                        DFN::Point_2D this_pnt{Intersection_AB[0]};
                        size_t IDX;
                        //if this intersection point is overlapping?
                        bool ft = this->Check_repetitive_pnt(this_pnt, IDX);

                        if (ft == true)
                        {
                            //not a new point, so IDX is the point NO
                            size_t idx_2;
                            bool pt = this->Insert_a_pnt_into_topo_line(Trace_topo_line_seg[i], this_pnt, idx_2);
                            if (pt == true)
                            {
                                Trace_topo_line_seg[i][idx_2].first = IDX;
                            }

                            size_t idx_3;
                            bool gy = this->Insert_a_pnt_into_topo_line(Trace_topo_line_seg[j], this_pnt, idx_3);
                            if (gy == true)
                            {
                                Trace_topo_line_seg[i][idx_3].first = IDX;
                            }
                        }
                        else
                        {
                            //must be a new point
                            size_t idx_2;
                            bool pt = this->Insert_a_pnt_into_topo_line(Trace_topo_line_seg[i], this_pnt, idx_2);
                            if (pt == true)
                            {
                                Pnt_sets.push_back(this_pnt);

                                Trace_topo_line_seg[i][idx_2].first = Pnt_sets.size() - 1; // a new point
                            }
                            else if (pt == false)
                            {
                                cout << "i = " << i << endl;
                                cout << "j = " << j << endl;
                                cout << "Intersection between traces should be a new point in this case 1!\n";
                                exit(0);
                            }

                            size_t idx_3;
                            bool gy = this->Insert_a_pnt_into_topo_line(Trace_topo_line_seg[j], this_pnt, idx_3);
                            if (gy == true)
                            {
                                Trace_topo_line_seg[j][idx_3].first = Pnt_sets.size() - 1; // a new point
                            }
                            else if (pt == false)
                            {
                                cout << "i = " << i << endl;
                                cout << "j = " << j << endl;
                                cout << "Intersection between traces should be a new point in this case 2!\n";
                                exit(0);
                            }
                        }
                    }
                    else if (cg == true && Intersection_AB.size() == 2)
                    {
                        //two traces overlapping
                        for (size_t k = 0; k < Intersection_AB.size(); ++k)
                        {
                            DFN::Point_2D this_pnt{Intersection_AB[k]};
                            size_t IDX;
                            //if this intersection point is overlapping?
                            bool ft = this->Check_repetitive_pnt(this_pnt, IDX);

                            if (ft == true)
                            {
                                //must be not a new point, so IDX is the point NO
                                size_t idx_2;
                                bool pt = this->Insert_a_pnt_into_topo_line(Trace_topo_line_seg[i], this_pnt, idx_2);
                                if (pt == true)
                                {
                                    Trace_topo_line_seg[i][idx_2].first = IDX;
                                }

                                size_t idx_3;
                                bool gy = this->Insert_a_pnt_into_topo_line(Trace_topo_line_seg[j], this_pnt, idx_3);
                                if (gy == true)
                                {
                                    Trace_topo_line_seg[j][idx_3].first = IDX;
                                }
                            }
                            else
                            {
                                cout << "if two traces overlap, the intersection points cannot be new points!\n";
                                exit(0);
                            }
                        }
                    }
                }
            }
        }

    for (size_t i = 0; i < Trace_topo_line_seg.size(); ++i)
        for (size_t pr = 0; pr < Trace_topo_line_seg[i].size() - 1; ++pr)
            if (this->Check_repetitive_line_seg(Trace_topo_line_seg[i][pr].first, Trace_topo_line_seg[i][pr + 1].first) == false)
                Topo_struct.push_back(make_pair(Trace_topo_line_seg[i][pr].first, Trace_topo_line_seg[i][pr + 1].first));
};

inline bool Polygon_convex_2D_with_traces::Insert_a_pnt_into_topo_line(std::vector<pair<size_t, DFN::Point_2D>> &topoline, DFN::Point_2D pnt, size_t &idx)
{
    for (size_t i = 0; i < topoline.size(); ++i)
    {
        bool y = pnt.If_overlap_with_another_pnt(topoline[i].second);
        if (y == true)
        {
            idx = topoline[i].first; // overlapping

            return false;
        }
    };

    for (size_t i = 1; i < topoline.size(); ++i)
    {
        double distancefront, distanceback;
        distancefront = (topoline[0].second.Coordinate - pnt.Coordinate).norm();
        distanceback = (topoline[0].second.Coordinate - topoline[i].second.Coordinate).norm();
        if (distancefront < distanceback)
        {
            idx = i;
            std::vector<pair<size_t, DFN::Point_2D>>::iterator it = topoline.begin();
            it = it + i;
            topoline.insert(it, 1, std::make_pair(1e10, pnt.Coordinate));

            return true;
        }
        else if (distancefront > distanceback && i == topoline.size() - 1)
        {
            idx = 1e10;
            return false;
        }
    }

    idx = 1e10;
    return false;
};

inline bool Polygon_convex_2D_with_traces::Check_repetitive_line_seg(const size_t A, const size_t B)
{
    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        if ((A == Topo_struct[i].first && B == Topo_struct[i].second) ||
            (B == Topo_struct[i].first && A == Topo_struct[i].second))
        {
            return true;
        }
    }
    return false;
};

inline bool Polygon_convex_2D_with_traces::Check_repetitive_pnt(const DFN::Point_2D G, size_t &idx)
{
    //std::vector<DFN::Point_2D> Pnt_sets;
    for (size_t i = 0; i < Pnt_sets.size(); ++i)
    {
        if (Pnt_sets[i].If_overlap_with_another_pnt(G) == true)
        {
            idx = i;
            return true;
        }
    }
    idx = 1e10;
    return false;
};

inline void Polygon_convex_2D_with_traces::Matlab_plot(string FileKey)
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

    for (size_t i = 0; i < Topo_struct.size(); ++i)
    {
        oss << "plot([Pnt_sets(" << Topo_struct[i].first + 1 << ",1), Pnt_sets(" << Topo_struct[i].second + 1 << ",1)],[Pnt_sets(" << Topo_struct[i].first + 1 << ",2), Pnt_sets(" << Topo_struct[i].second + 1 << ",2)], 'LineWidth',2);\n";
        oss << "hold on;\n";
    }

    oss.close();
};

}; // namespace DFN