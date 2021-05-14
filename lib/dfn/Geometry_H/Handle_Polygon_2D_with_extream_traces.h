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

////////////////// not robust algorithm

namespace DFN
{
class Handle_Polygon_2D_with_extream_traces
{
public:
    std::vector<Vector3d> Verts_handled;

public:
    Handle_Polygon_2D_with_extream_traces(std::vector<Vector3d> Verts, vector<DFN::Line_seg_2D> Traces, const double subseg);
    void Insert_a_point(std::vector<Vector3d> &Verts_sort_sub, const Vector3d P_2D, const double distance_foot);
    void Matlab_plot(string FileKey, vector<DFN::Line_seg_2D> Traces);
};

//------------------------------

inline Handle_Polygon_2D_with_extream_traces::Handle_Polygon_2D_with_extream_traces(std::vector<Vector3d> Verts, vector<DFN::Line_seg_2D> Traces, const double subseg)
{
    std::vector<std::vector<Vector3d>> Verts_sort(Verts.size());
    for (size_t i = 0; i < Verts.size(); ++i)
    {
        Verts_sort[i].resize(2);
        Verts_sort[i][0] = Verts[i];
        Verts_sort[i][1] = Verts[(i + 1) % Verts.size()];
    }

    for (size_t i = 0; i < Traces.size(); ++i)
    {
        if (Traces[i].If_is_a_point() == true) // if trace is a point
        {
            DFN::Point_2D ThisPnt{Traces[i].Point[0]};
            //knowing that this point must be within or lie on the polygon
            for (size_t j = 0; j < Verts.size(); ++j)
            {
                DFN::Distance Disy{};
                DFN::Line_seg_2D Line_seg1{Vector2d{Verts[j](0), Verts[j](1)}, Vector2d{Verts[(j + 1) % Verts.size()](0), Verts[(j + 1) % Verts.size()](1)}};
                bool lieson;
                double distance_k = Disy.Distance_between_point_lineseg(ThisPnt, Line_seg1, lieson);

                if (lieson == true)
                {
                    //do nothing
                }
                else
                {
                    if (distance_k < subseg)
                    {
                        // need to know the normal foot
                        Vector2d foot = ThisPnt.Perpend_foot_on_a_line_seg(Line_seg1.Point);
                        // now insert a point to the polygon frame
                        this->Insert_a_point(Verts_sort[j], Vector3d{foot[0], foot[1], 0}, distance_k);
                    }
                    else
                    {
                        //do noting
                    }
                }
            }
        }
        else // if trace is not a point
        {
            //a trace has two ends
            for (size_t k = 0; k < 2; ++k)
            {
                DFN::Point_2D ThisPnt{Traces[i].Point[k]};

                for (size_t j = 0; j < Verts.size(); ++j)
                {
                    DFN::Distance Disy{};
                    DFN::Line_seg_2D Line_seg1{Vector2d{Verts[j](0), Verts[j](1)}, Vector2d{Verts[(j + 1) % Verts.size()](0), Verts[(j + 1) % Verts.size()](1)}};
                    bool lieson;
                    double distance_k = Disy.Distance_between_point_lineseg(ThisPnt, Line_seg1, lieson);

                    if (lieson == true)
                    {
                        //do nothing
                    }
                    else
                    {
                        if (distance_k < subseg)
                        {
                            // need to know the normal foot
                            Vector2d foot = ThisPnt.Perpend_foot_on_a_line_seg(Line_seg1.Point);
                            // now insert a point to the polygon frame
                            this->Insert_a_point(Verts_sort[j], Vector3d{foot[0], foot[1], 0}, distance_k);
                        }
                        else
                        {
                            //do noting
                        }
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < Verts_sort.size(); ++i)
    {
        for (size_t j = 0; j < Verts_sort[i].size() - 1; ++j)
        {
            Verts_handled.push_back(Verts_sort[i][j]);
            //cout << Verts_handled.size() << endl;
        }
    }
};

void Handle_Polygon_2D_with_extream_traces::Insert_a_point(std::vector<Vector3d> &Verts_sort_sub, const Vector3d P_2D, const double distance_foot)
{
    std::vector<Vector3d> Verts_sort_sub_2(1);
    Verts_sort_sub_2[0] << Verts_sort_sub[0](0), Verts_sort_sub[0](1), 0;

    double distance_h = (Vector3d{P_2D[0], P_2D[1], 0} - Verts_sort_sub_2[0]).norm();

    for (size_t i = 1; i < Verts_sort_sub.size(); ++i)
    {
        double distance_r = (Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0} - Verts_sort_sub_2[0]).norm();

        if (abs(distance_r) < abs(distance_h))
        {
            Verts_sort_sub_2.push_back(Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0});

            if (i == Verts_sort_sub.size() - 1)
            {
                cout << "This point is not above or below the line seg!\nIn class: Handle_Polygon_2D_with_extream_traces, function: Insert_a_point\n";
                exit(0);
            }
        }
        else
        {
            double distance_front = abs((Vector3d{P_2D[0], P_2D[1], 0} - Vector3d{Verts_sort_sub[i - 1][0], Verts_sort_sub[i - 1][1], 0}).norm());
            double distance_back = abs((Vector3d{P_2D[0], P_2D[1], 0} - Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0}).norm());
            //cout << 1 << endl;
            if (distance_front < 0.1 || distance_back < 0.1)
            {
                Verts_sort_sub_2.push_back(Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0});
            }
            else
            {
                if (distance_front <= 0.5 * distance_foot && distance_back <= 0.5 * distance_foot)
                {
                    Verts_sort_sub_2.push_back(Vector3d{P_2D[0], P_2D[1], 0});
                    Verts_sort_sub_2.push_back(Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0});
                }
                else if (distance_front > 0.5 * distance_foot && distance_back <= 0.5 * distance_foot)
                {
                    //push points
                    Line_seg_2D Line_f{Vector2d{Verts_sort_sub[i - 1][0], Verts_sort_sub[i - 1][1]}, Vector2d{P_2D[0], P_2D[1]}};
                    std::vector<Vector2d> KKY = Line_f.Exact_pnts_along_this_line_seg_evenly(0.5 * distance_foot);

                    Verts_sort_sub_2.push_back(Vector3d{KKY[KKY.size() - 1](0), KKY[KKY.size() - 1](1), 0});

                    Verts_sort_sub_2.push_back(Vector3d{P_2D[0], P_2D[1], 0});
                    Verts_sort_sub_2.push_back(Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0});
                }
                else if (distance_front <= 0.5 * distance_foot && distance_back > 0.5 * distance_foot)
                {
                    Verts_sort_sub_2.push_back(Vector3d{P_2D[0], P_2D[1], 0});
                    //push points
                    Line_seg_2D Line_f{Vector2d{P_2D[0], P_2D[1]}, Vector2d{Verts_sort_sub[i][0], Verts_sort_sub[i][1]}};
                    std::vector<Vector2d> KKY = Line_f.Exact_pnts_along_this_line_seg_evenly(0.5 * distance_foot);

                    Verts_sort_sub_2.push_back(Vector3d{KKY[0](0), KKY[0](1), 0});

                    Verts_sort_sub_2.push_back(Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0});
                }
                else
                {

                    Line_seg_2D Line_f{Vector2d{Verts_sort_sub[i - 1][0], Verts_sort_sub[i - 1][1]}, Vector2d{P_2D[0], P_2D[1]}};
                    std::vector<Vector2d> KKY = Line_f.Exact_pnts_along_this_line_seg_evenly(0.5 * distance_foot);

                    Verts_sort_sub_2.push_back(Vector3d{KKY[KKY.size() - 1](0), KKY[KKY.size() - 1](1), 0});

                    Verts_sort_sub_2.push_back(Vector3d{P_2D[0], P_2D[1], 0});
                    //----------------

                    Line_seg_2D Line_ts{Vector2d{P_2D[0], P_2D[1]}, Vector2d{Verts_sort_sub[i][0], Verts_sort_sub[i][1]}};
                    std::vector<Vector2d> KKe = Line_ts.Exact_pnts_along_this_line_seg_evenly(0.5 * distance_foot);

                    Verts_sort_sub_2.push_back(Vector3d{KKe[0](0), KKe[0](1), 0});

                    Verts_sort_sub_2.push_back(Vector3d{Verts_sort_sub[i][0], Verts_sort_sub[i][1], 0});
                }
            }
        }
    }

    Verts_sort_sub = Verts_sort_sub_2;
};

void Handle_Polygon_2D_with_extream_traces::Matlab_plot(string FileKey, vector<DFN::Line_seg_2D> Traces)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";

    for (size_t i = 0; i < Verts_handled.size(); ++i)
    {
        oss << "plot([";
        oss << Verts_handled[i](0) << ", " << Verts_handled[(i + 1) % Verts_handled.size()](0) << "], [";
        oss << Verts_handled[i](1) << ", " << Verts_handled[(i + 1) % Verts_handled.size()](1) << "], '-*', 'linewidth', 1);\n";
        oss << "hold on;\n";
    }

    for (size_t i = 0; i < Traces.size(); ++i)
    {
        oss << "plot([";
        oss << Traces[i].Point[0](0) << ", " << Traces[i].Point[1](0) << "], [";
        oss << Traces[i].Point[0](1) << ", " << Traces[i].Point[1](1) << "], '-o', 'linewidth', 2);\n";
        oss << "hold on;\n";
    }

    oss.close();
};
}; // namespace DFN