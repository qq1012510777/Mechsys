#pragma once
#include <algorithm>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Fade_2D
#include "DFN_fade2d_WL.h"

namespace DFN
{

class DFN_mesh
{
public:
    std::vector<std::pair<Vector4d, Vector6d>> Rota_angle;
    ///< 1. the rotation angle and the center of the
    // fracture 
    //< 2.
    // the normal of fracture and the
    // vector perpendicular to the normal (in xy plane)

    std::vector<std::vector<Vector3d>> JXY;
    ///< JXY.size() is the number of fractures,
    //JXY[0].size() is the number of vertices of
    //1st fractures (also with ID NO 1)

    std::vector<std::vector<Vector3d>> JXY_3D;
    ///<

    std::vector<size_t> JXY_frac_ID;
    ///<

    std::vector<std::vector<Vector6d>> JM;
    ///< each Vector3d is the ID values of
    //the three nodes

    std::vector<std::map<size_t, double>> JB_1;
    ///< 1st kind of BC

    std::vector<std::map<std::pair<size_t, size_t>, Vector5d>> JB_2;
    ///< 2nd kind of BC: ele no, edge no (local), cos(theta_x), cos(theta_y), q1, q2, q3

    std::vector<std::vector<std::vector<size_t>>> Trace_Node_overall;
    ///< node NOs of each traces in each fracture

    std::vector<std::vector<Vector3d>> Coe_Matr_guide;
    ///< (0): an index referring if this is
    //a repetitive point;
    //(1): i; (2): j;(if repetitive)

    std::vector<std::vector<Vector3d>> Coe_Matr_guide_m;
    ///< (0): an index referring if this is
    //a repetitive point;
    //(1): i; (2): j;(if repetitive)

    size_t NO_all_pnts;
    ///<

    size_t NO_elements;
    ///<

    size_t NO_Nodes_p;

    std::vector<size_t> NO_Nodes_p_each_frac;

    std::vector<std::vector<int>> Listofclusters_mesh_only;
    ///<

    std::vector<std::map<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>>> Inlet;
    ///< the ith frac, and the jth element,
    //< first vector6d: 2d pnt 1 and 2;
    ///< additionally, Inlet.second.first(2) = edge no
    //< second vector6d: 3d pnt 1 and 2

    std::vector<std::map<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>>> Outlet;
    ///<

    std::vector<std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>>> neigh_shared;
    ///< i = i fracture, j = j element, size_t = neighbor ele ID,, Vector2d = shared edge no of j ele, shared edge no of neighbor ele

    //std::vector<std::vector<Vector6d>> JB_2;
    ///< 2nd kind of BC

    // function
public:
    DFN_mesh(const DFN::Domain dom,
             const double avg_ele_len,
             const double ratio_H,
             string percolation_direction);
    //constructor

    void Find_pnt_inside_region(const DFN::Domain dom,
                                const std::vector<Vector3d> tem_verts_trim,
                                const size_t i,
                                const size_t P_idx,
                                const size_t he,
                                double &x_ik,
                                double &y_ik);
    // <

    void Rotation_3D_frac_to_2D(DFN::Domain dom,
                                std::vector<Vector3d> &tem_verts_trim_AAA,
                                std::map<std::pair<double, double>, int> &MapPnt,
                                std::vector<Vector2d> &SegIDtoTD,
                                Vector3d &temp3,
                                Vector3d &Normal_frac,
                                const size_t i,
                                const size_t P_idx,
                                const size_t he,
                                const double d_In,
                                std::vector<std::pair<Vector3d, Vector3d>> &Connection_traces);
    //<

    void Matlab_Plot_DFN_mesh(string FileKey,
                              const DFN::Domain dom,
                              const std::vector<std::vector<Vector3d>> JXY,
                              const std::vector<std::vector<Vector6d>> JM,
                              std::vector<std::map<size_t, double>> JB_1,
                              std::vector<std::map<std::pair<size_t, size_t>, Vector5d>> JB_2,
                              const std::vector<std::vector<std::vector<size_t>>> Trace_Node_overall,
                              const size_t P_idx);
    //<

    void CXX_Triangle_code(string FileKey,
                           const size_t i,
                           const std::vector<Vector3d> tem_verts_trim,
                           const std::vector<std::vector<Vector3d>> Seg_pnt,
                           const double x_ik,
                           const double y_ik,
                           const double element_max);

    void CXX_Fade2D_code(string FileKey,
                         const size_t i,
                         std::map<std::pair<double, double>, int> MapPnt,
                         std::vector<Vector2d> SegIDtoID,
                         const std::vector<Vector3d> tem_verts_trim_AAA);

    void Matlab_plot_2D_frac(string FileKey,
                             const size_t i,
                             const std::vector<Vector3d> tem_verts_trim,
                             const std::vector<std::vector<Vector3d>> Seg_pnt);

    void Address_Repetitive_Node(std::vector<std::vector<Vector3d>> &Coe_Matr_guide,
                                 std::vector<std::vector<Vector3d>> &Coe_Matr_guide_m,
                                 const std::vector<std::vector<Vector3d>> JXY,
                                 size_t &overall_matrix_dimension);

    void Remove_unnecessary_fractures(size_t Cluster_ID,
                                      DFN::Domain dom);
};

///************************************************************///
inline DFN_mesh::DFN_mesh(const DFN::Domain dom,
                          const double avg_ele_len,
                          const double ratio_H,
                          string percolation_direction)
{
    NO_Nodes_p = 0;
    Listofclusters_mesh_only.resize(dom.Listofclusters.size());
    for (size_t i = 0; i < dom.Listofclusters.size(); ++i)
    {
        Listofclusters_mesh_only[i].resize(dom.Listofclusters[i].size());
        for (size_t j = 0; j < Listofclusters_mesh_only[i].size(); ++j)
        {
            Listofclusters_mesh_only[i][j] = dom.Listofclusters[i][j];
        }
    }

    NO_elements = 0;
    if (JM.size() >= 1 || JXY.size() >= 1 || JXY_3D.size() >= 1)
    {
        std::cout << "Error in Class 'Mesh_frac', size of JXY, JXY_3D or JM should be initialized!\n";
    }

    size_t P_idx;
    if (percolation_direction == "x")
        P_idx = 0;
    else if (percolation_direction == "y")
        P_idx = 1;
    else if (percolation_direction == "z")
        P_idx = 2;
    else
    {
        std::cout << "Error! Cannot recognize the per-defined percolation direction!\n";
        exit(0);
    };

    //std::vector<std::vector<size_t>> Listofclusters;

    for (size_t he = 0; he < dom.Percolation_cluster[P_idx].size(); ++he)
    {
        //std::cout << "000\n";
        //dom.Listofclusters[dom.Percolation_cluster[P_idx][he]], this is the percolation cluster
        Remove_unnecessary_fractures(dom.Percolation_cluster[P_idx][he], dom);
        //std::cout << "111\n";
        for (size_t i = 0; i < dom.Listofclusters[dom.Percolation_cluster[P_idx][he]].size(); ++i)
        {
            //std::cout << "111\n";
            if (Listofclusters_mesh_only[dom.Percolation_cluster[P_idx][he]][i] != -1)
            {
                int frac_ID_k = Listofclusters_mesh_only[dom.Percolation_cluster[P_idx][he]][i];
                JXY_frac_ID.push_back((size_t)(frac_ID_k));
                std::vector<Vector3d> JXY_1;
                std::vector<Vector3d> JXY_3D_1;
                std::vector<Vector6d> JM_1;
                std::pair<Vector4d, Vector6d> Rota_angle_1;
                std::vector<std::vector<size_t>> Trace_Node_A;
                std::map<size_t, double> JBA_1;                      ///< 1st kind of BC
                std::map<std::pair<size_t, size_t>, Vector5d> JBA_2; ///< 2nd kind of BC
                std::vector<std::pair<Vector3d, Vector3d>> Connection_traces;
                //std::vector<Vector6d> JBA_2; ///< 2nd kind of BC
                //std::vector<Vector3d> JXY_inserted_pnt_1; // this variate can be removed

                double d_In = avg_ele_len; ///< edge length of equilateral triangle (the ideal triangular element)

                // rotation from 3D to 2D starts//
                Vector3d Normal_frac; ///< the normal vector of this fracture
                Vector3d temp3;       ///< a vector perpendicular to the normal vector of fracture, but this one is lying in the x-y plane
                std::map<std::pair<double, double>, int> MapPnt;
                std::vector<Vector2d> SegIDtoID;
                std::vector<Vector3d> tem_verts_trim_AAA;
                Rotation_3D_frac_to_2D(dom, tem_verts_trim_AAA, MapPnt, SegIDtoID, temp3, Normal_frac, i, P_idx, he, d_In, Connection_traces); ///< center is (0, 0)
                //std::cout << "111;\n";
                // rotation from 3D to 2D ends//
                //CXX_Fade2D_code("tdfn_FadeCXX", i, MapPnt, SegIDtoID, tem_verts_trim_AAA);

                //ratio_H; // the ratio of max element edge length to a sub_segment of frac boundary or intersection trace
                double max_len_ele = avg_ele_len * ratio_H;

                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> neigh_shared_A;

                size_t NO_Nodes_p_A = 0;
                DFN::FADE2D FadeMesh(MapPnt, SegIDtoID, 29, 0.5, max_len_ele, neigh_shared_A, NO_Nodes_p_A);

                //store 2D mesh vertexes in JXY_1, JXY_3D_1
                JXY_1.resize(FadeMesh.JXY.size());
                JXY_3D_1.resize(FadeMesh.JXY.size());
                for (size_t j = 0; j < JXY_1.size(); ++j)
                {
                    Vector3d temyh;
                    temyh << FadeMesh.JXY[j](0), FadeMesh.JXY[j](1), 0;
                    JXY_1[j] = temyh;
                    JXY_3D_1[j] = temyh;
                }

                Rota_angle_1.first(0) = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Dip_angle * M_PI / 180;
                Rota_angle_1.first(1) = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center(0);
                Rota_angle_1.first(2) = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center(1);
                Rota_angle_1.first(3) = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center(2);

                Rota_angle_1.second(0) = Normal_frac(0);
                Rota_angle_1.second(1) = Normal_frac(1);
                Rota_angle_1.second(2) = Normal_frac(2);
                Rota_angle_1.second(3) = temp3(0);
                Rota_angle_1.second(4) = temp3(1);
                Rota_angle_1.second(5) = temp3(2);

                if (abs(Normal_frac(0)) < 0.000001 && abs(Normal_frac(1)) < 0.000001 && abs(Normal_frac(2)) < 0.000001)
                {
                    for (size_t j = 0; j < JXY_1.size(); ++j)
                    {
                        JXY_3D_1[j] += dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center;
                    }
                }
                else
                {
                    double R_angle_temp1 = 0;
                    double x_temp = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Dip_angle; ///it is better to create a new variable to represent the dip angle, because debuging shows direct use of 'Dip_angle' to calculate rotation angle leads wrong output
                    R_angle_temp1 = x_temp * M_PI / 180;

                    Quaternion_t Q_axis_1;
                    NormalizeRotation(R_angle_temp1, temp3, Q_axis_1);

                    for (size_t j = 0; j < JXY_1.size(); ++j)
                    {
                        Vector3d temp4;
                        Vector3d temp5;
                        temp5 = JXY_3D_1[j];
                        Rotation(temp5, Q_axis_1, temp4);
                        temp5 = temp4 + dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center;
                        JXY_3D_1[j] = temp5;
                    };
                }

                JM_1.resize(FadeMesh.JM.size());
                NO_elements += FadeMesh.JM.size();
                //int kjj = 0;
                for (size_t j = 0; j < FadeMesh.JM.size(); ++j)
                {
                    Vector6d temyjk;
                    temyjk << FadeMesh.JM[j](0, 0),
                        FadeMesh.JM[j](0, 1),
                        FadeMesh.JM[j](0, 2),
                        FadeMesh.JM[j](0, 3),
                        FadeMesh.JM[j](0, 4),
                        FadeMesh.JM[j](0, 5);
                    /*
                        std::vector<Vector3d> A(3);
                        A[0] = JXY_1[mesh.Cells[j]->V[0]->ID];
                        A[1] = JXY_1[mesh.Cells[j]->V[1]->ID];
                        A[2] = JXY_1[mesh.Cells[j]->V[2]->ID]; 
                        if (Area_of_a_convex(A) > element_max && kjj < 1)
                        {
                            std::cout << "The Fracture ID: " << i << "! Found large element!!!\n";
                            kjj++;
                        }
                        */
                    JM_1[j] = temyjk;
                }

                ///< Vector6d Model_domain;
                ///< 0 Top-zmax, 1 bottom-zmin, 2 front-ymin, 3 back-ymax, 4 left-xmin, 5 right-xmax
                ///< identify domain and element boundary nodes
                std::vector<double> BC_DOMAIN_FLAG = {1, 1, 2, 2, 2, 2}; ///< indicate top and bottom domain faces are of 1st BCs
                std::vector<double> BC_DOMAIN_VALUE_1st = {100, 20, 0, 0, 0, 0};
                std::vector<double> BC_DOMAIN_VALUE_2nd = {0, 0, 0, 0, 0, 0}; // impermeable

                //first, domain boundaries are of 1st kind of BC

                for (size_t j = 0; j < NO_Nodes_p_A; ++j)
                {
                    Vector3d PNT_3D;
                    PNT_3D = JXY_3D_1[j];
                    std::vector<bool> B_Domain(6);
                    size_t vert_ID = j;
                    double h_tmp;

                    if (BC_DOMAIN_FLAG[0] == 1)
                    {
                        Vector3d A_domain_top, B_domain_top, C_domain_top, D_domain_top;
                        A_domain_top << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                        B_domain_top << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                        C_domain_top << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                        D_domain_top << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                        B_Domain[0] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_top, B_domain_top, C_domain_top, D_domain_top, PNT_3D);
                        if (B_Domain[0] == true)
                            h_tmp = BC_DOMAIN_VALUE_1st[0];
                    }

                    if (BC_DOMAIN_FLAG[1] == 1)
                    {
                        Vector3d A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom;
                        A_domain_bottom << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                        B_domain_bottom << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                        C_domain_bottom << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                        D_domain_bottom << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                        B_Domain[1] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom, PNT_3D);
                        if (B_Domain[1] == true)
                            h_tmp = BC_DOMAIN_VALUE_1st[1];
                    }

                    if (BC_DOMAIN_FLAG[2] == 1)
                    {
                        Vector3d A_domain_front, B_domain_front, C_domain_front, D_domain_front;
                        A_domain_front << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                        B_domain_front << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                        C_domain_front << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                        D_domain_front << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                        B_Domain[2] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_front, B_domain_front, C_domain_front, D_domain_front, PNT_3D);
                        if (B_Domain[2] == true)
                            h_tmp = BC_DOMAIN_VALUE_1st[2];
                    }

                    if (BC_DOMAIN_FLAG[3] == 1)
                    {
                        Vector3d A_domain_back, B_domain_back, C_domain_back, D_domain_back;
                        A_domain_back << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                        B_domain_back << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                        C_domain_back << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                        D_domain_back << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                        B_Domain[3] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_back, B_domain_back, C_domain_back, D_domain_back, PNT_3D);
                        if (B_Domain[3] == true)
                            h_tmp = BC_DOMAIN_VALUE_1st[3];
                    }

                    if (BC_DOMAIN_FLAG[4] == 1)
                    {
                        Vector3d A_domain_left, B_domain_left, C_domain_left, D_domain_left;
                        A_domain_left << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                        B_domain_left << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                        C_domain_left << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                        D_domain_left << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                        B_Domain[4] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_left, B_domain_left, C_domain_left, D_domain_left, PNT_3D);
                        if (B_Domain[4] == true)
                            h_tmp = BC_DOMAIN_VALUE_1st[4];
                    }

                    if (BC_DOMAIN_FLAG[5] == 1)
                    {
                        Vector3d A_domain_right, B_domain_right, C_domain_right, D_domain_right;
                        A_domain_right << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                        B_domain_right << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                        C_domain_right << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                        D_domain_right << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                        B_Domain[5] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_right, B_domain_right, C_domain_right, D_domain_right, PNT_3D);
                        if (B_Domain[5] == true)
                            h_tmp = BC_DOMAIN_VALUE_1st[5];
                    }

                    if (B_Domain[0] == true ||
                        B_Domain[1] == true ||
                        B_Domain[2] == true ||
                        B_Domain[3] == true ||
                        B_Domain[4] == true ||
                        B_Domain[5] == true)
                    {
                        JBA_1.insert(std::pair<size_t, double>(vert_ID, h_tmp));
                    };
                }

                //second, domain boundaries are of 2nd kind of BC
                for (size_t j = 0; j < JM_1.size(); ++j)
                {
                    std::vector<std::map<size_t, size_t>> MapBC2nd(6); // 6 means six faces
                    std::vector<int> MapSize = {0, 0, 0, 0, 0, 0};
                    for (size_t k = 0; k < 6; ++k)
                    {
                        Vector3d PNT_3D;
                        PNT_3D = JXY_3D_1[JM_1[j][k]];
                        std::vector<bool> B_Domain(6);
                        //size_t vert_ID = JM_1[j][k];
                        //double h_tmp;

                        if (BC_DOMAIN_FLAG[0] == 2)
                        {
                            Vector3d A_domain_top, B_domain_top, C_domain_top, D_domain_top;
                            A_domain_top << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                            B_domain_top << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                            C_domain_top << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                            D_domain_top << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                            bool po = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_top, B_domain_top, C_domain_top, D_domain_top, PNT_3D);
                            if (po == true)
                            {
                                std::pair<std::map<size_t, size_t>::iterator, bool> ret = MapBC2nd[0].insert(std::pair<size_t, size_t>(k, 0));
                                if (ret.second)
                                {
                                    MapSize[0]++;
                                }
                            }
                        }

                        if (BC_DOMAIN_FLAG[1] == 2)
                        {
                            Vector3d A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom;
                            A_domain_bottom << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                            B_domain_bottom << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                            C_domain_bottom << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                            D_domain_bottom << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                            bool po = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom, PNT_3D);
                            if (po == true)
                            {
                                std::pair<std::map<size_t, size_t>::iterator, bool> ret = MapBC2nd[1].insert(std::pair<size_t, size_t>(k, 1));
                                if (ret.second)
                                {
                                    MapSize[1]++;
                                }
                            }
                        }

                        if (BC_DOMAIN_FLAG[2] == 2)
                        {
                            Vector3d A_domain_front, B_domain_front, C_domain_front, D_domain_front;
                            A_domain_front << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                            B_domain_front << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                            C_domain_front << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                            D_domain_front << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                            bool po = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_front, B_domain_front, C_domain_front, D_domain_front, PNT_3D);
                            if (po == true)
                            {
                                std::pair<std::map<size_t, size_t>::iterator, bool> ret = MapBC2nd[2].insert(std::pair<size_t, size_t>(k, 2));
                                if (ret.second)
                                {
                                    MapSize[2]++;
                                }
                            }
                        }

                        if (BC_DOMAIN_FLAG[3] == 2)
                        {
                            Vector3d A_domain_back, B_domain_back, C_domain_back, D_domain_back;
                            A_domain_back << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                            B_domain_back << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                            C_domain_back << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                            D_domain_back << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                            bool po = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_back, B_domain_back, C_domain_back, D_domain_back, PNT_3D);
                            if (po == true)
                            {
                                std::pair<std::map<size_t, size_t>::iterator, bool> ret = MapBC2nd[3].insert(std::pair<size_t, size_t>(k, 3));
                                if (ret.second)
                                {
                                    MapSize[3]++;
                                }
                            }
                        }

                        if (BC_DOMAIN_FLAG[4] == 2)
                        {
                            Vector3d A_domain_left, B_domain_left, C_domain_left, D_domain_left;
                            A_domain_left << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                            B_domain_left << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                            C_domain_left << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                            D_domain_left << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                            bool po = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_left, B_domain_left, C_domain_left, D_domain_left, PNT_3D);
                            if (po == true)
                            {
                                std::pair<std::map<size_t, size_t>::iterator, bool> ret = MapBC2nd[4].insert(std::pair<size_t, size_t>(k, 4));
                                if (ret.second)
                                {
                                    MapSize[4]++;
                                }
                            }
                        }

                        if (BC_DOMAIN_FLAG[5] == 2)
                        {
                            Vector3d A_domain_right, B_domain_right, C_domain_right, D_domain_right;
                            A_domain_right << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                            B_domain_right << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                            C_domain_right << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                            D_domain_right << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                            bool po = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_right, B_domain_right, C_domain_right, D_domain_right, PNT_3D);
                            if (po == true)
                            {
                                std::pair<std::map<size_t, size_t>::iterator, bool> ret = MapBC2nd[5].insert(std::pair<size_t, size_t>(k, 5));
                                if (ret.second)
                                {
                                    MapSize[5]++;
                                }
                            }
                        };
                    }
                    //std::vector<std::map<size_t, size_t>> MapBC2nd(6); // 6 means six faces
                    //std::cout << "side boundary:\n";
                    for (size_t k = 0; k < MapBC2nd.size(); ++k)
                    {
                        if (BC_DOMAIN_FLAG[k] == 2)
                        {
                            if (MapSize[k] == 3)
                            {
                                //push std::map<std::pair<size_t, size_t>, Vector4d> JBA_2; ///< 2nd kind of BC
                                std::map<size_t, size_t>::iterator it_A, itPLus_1_A, itPLus_2_A, itEnd_A;
                                it_A = MapBC2nd[k].begin();
                                itPLus_1_A = MapBC2nd[k].begin();
                                itPLus_2_A = MapBC2nd[k].begin();
                                itPLus_1_A++;
                                itPLus_2_A++;
                                itPLus_2_A++;
                                itEnd_A = MapBC2nd[k].end();

                                size_t EDGE_NO = 0;
                                /*
                                        std::cout << "it_A->first: " << it_A->first << "; "
                                        << "itPLus_1_A->first: " << itPLus_1_A->first << "\n";
                                    */
                                if ((it_A->first == 0) && (itPLus_1_A->first == 1) && (itPLus_2_A->first == 2))
                                {
                                    EDGE_NO = 0;
                                }
                                else if ((it_A->first == 2) && (itPLus_1_A->first == 3) && (itPLus_2_A->first == 4))
                                {
                                    EDGE_NO = 1;
                                }
                                else if ((it_A->first == 0) && (itPLus_1_A->first == 4) && (itPLus_2_A->first == 5))
                                {
                                    EDGE_NO = 2;
                                }
                                else
                                {
                                    std::cout << "Cannot define Edge when identifying boundary vertices!\n";
                                    exit(0);
                                }
                                std::pair<size_t, size_t> tmp_pair = std::make_pair(j, EDGE_NO); // element NO, edge NO
                                Vector5d tmp_yu;
                                double cos_theta_x = 0, cos_theta_y = 0;

                                tmp_yu << cos_theta_x, cos_theta_y, 0, 0, 0;

                                JBA_2.insert(std::pair<std::pair<size_t, size_t>, Vector5d>(tmp_pair, tmp_yu));
                                //std::pair<std::map<std::pair<size_t, size_t>, Vector5d>::iterator, bool> ret = JBA_2.insert(std::pair<std::pair<size_t, size_t>, Vector5d>(tmp_pair, tmp_yu));
                                /*
                                    if (ret.second)
                                    {
                                        std::cout << "the " << k << "th face, edge->" << EDGE_NO << ", element NO->" << j<< " : ";
                                        while (it_A != itEnd_A)
                                        {
                                            std::cout << it_A->first << ", [";
                                            std::cout << JXY_3D_1[JM_1[j](it_A->first)] << "], ";
                                            it_A++;
                                        }
                                        std::cout << "\n\n\n";
                                    }*/
                            }
                        }
                    };
                }

                ///< third, fracture boundaries are of 2nd kind of BC, especially q = 0
                ///< so if fracture boundaries are also belonging to domain boundaries, then they should not be counted again
                //std::vector<Vector3d> tem_verts_trim_AAA;
                for (size_t j = 0; j < JM_1.size(); ++j)
                {
                    std::vector<std::map<size_t, size_t>> MapBC2nd(tem_verts_trim_AAA.size()); // how many edges enclosing a fracture
                    std::vector<int> MapSize(tem_verts_trim_AAA.size());
                    for (size_t k = 0; k < MapSize.size(); ++k)
                        MapSize[k] = 0;

                    for (size_t k = 0; k < 6; ++k)
                    {
                        Vector3d PNT_2D;
                        size_t vert_NO = JM_1[j][k];
                        PNT_2D = JXY_1[vert_NO];
                        for (size_t h0 = 0; h0 < tem_verts_trim_AAA.size(); ++h0)
                        {
                            bool Tf = Determine_IF_a_pnt_lies_in_a_2D_line_segment(tem_verts_trim_AAA[h0], tem_verts_trim_AAA[(h0 + 1) % (tem_verts_trim_AAA.size())], PNT_2D);
                            if (Tf == true)
                            {
                                std::pair<std::map<size_t, size_t>::iterator, bool> ret = MapBC2nd[h0].insert(std::pair<size_t, size_t>(k, h0));
                                if (ret.second)
                                {
                                    MapSize[h0]++;
                                }
                            }
                        };
                    };
                    //std::cout << "frac boundary:\n";
                    for (size_t k = 0; k < MapBC2nd.size(); ++k)
                    {
                        if (MapSize[k] == 3)
                        {
                            //push std::map<std::pair<size_t, size_t>, Vector4d> JBA_2; ///< 2nd kind of BC
                            std::map<size_t, size_t>::iterator it_A, itPLus_1_A, itPLus_2_A, itEnd_A;
                            it_A = MapBC2nd[k].begin();
                            itPLus_1_A = MapBC2nd[k].begin();
                            itPLus_2_A = MapBC2nd[k].begin();
                            itPLus_1_A++;
                            itPLus_2_A++;
                            itPLus_2_A++;
                            itEnd_A = MapBC2nd[k].end();

                            size_t EDGE_NO = 0;
                            /*
                                std::cout << "it_A->first: " << it_A->first << "; "
                                        << "itPLus_1_A->first: " << itPLus_1_A->first << "\n";
                                */
                            if ((it_A->first == 0) && (itPLus_1_A->first == 1) && (itPLus_2_A->first == 2))
                            {
                                EDGE_NO = 0;
                            }
                            else if ((it_A->first == 2) && (itPLus_1_A->first == 3) && (itPLus_2_A->first == 4))
                            {
                                EDGE_NO = 1;
                            }
                            else if ((it_A->first == 0) && (itPLus_1_A->first == 4) && (itPLus_2_A->first == 5))
                            {
                                EDGE_NO = 2;
                            }
                            else
                            {
                                std::cout << "Cannot define Edge when identifying boundary vertices!\n";
                                exit(0);
                            }
                            std::pair<size_t, size_t> tmp_pair = std::make_pair(j, EDGE_NO); // element NO, edge NO

                            //--also need to make sure the the edge is not adjacent to the domain boundary of 1st boundary condition
                            size_t Local_node_0 = 0, Local_node_2 = 0;
                            if (EDGE_NO == 0)
                            {
                                Local_node_0 = 0;
                                Local_node_2 = 2;
                            }
                            else if (EDGE_NO == 1)
                            {
                                Local_node_0 = 2;
                                Local_node_2 = 4;
                            }
                            else if (EDGE_NO == 2)
                            {
                                Local_node_0 = 4;
                                Local_node_2 = 0;
                            }
                            else
                            {
                                std::cout << "Error! Undefined Edge NO!\n";
                                exit(0);
                            }
                            Vector3d tmp_A1, tmp_B1;
                            tmp_A1 = JXY_3D_1[JM_1[j][Local_node_0]];
                            tmp_B1 = JXY_3D_1[JM_1[j][Local_node_2]];

                            bool itt[2];
                            itt[0] = false;
                            itt[1] = false;
                            //std::cout << "ELeNO: " << j << ", Edge No: " << EDGE_NO <<"; Now, checking if it is belonging to top and bottom boundaries!\n";
                            //std::cout << "The two pnts: " << tmp_A1 << ", " << tmp_B1 << std::endl;
                            for (size_t jk = 0; jk < 1; ++jk)
                            {
                                if (BC_DOMAIN_FLAG[0] == 1)
                                {
                                    Vector3d A_domain_top, B_domain_top, C_domain_top, D_domain_top;
                                    A_domain_top << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                                    B_domain_top << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                                    C_domain_top << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                                    D_domain_top << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                                    itt[0] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_top, B_domain_top, C_domain_top, D_domain_top, tmp_A1);
                                    itt[1] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_top, B_domain_top, C_domain_top, D_domain_top, tmp_B1);
                                    //std::cout << "TOP: " << itt[0] <<", " << itt[1] << std::endl;
                                    if (itt[0] == true && itt[1] == true)
                                    {
                                        break;
                                    };
                                }

                                if (BC_DOMAIN_FLAG[1] == 1)
                                {
                                    Vector3d A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom;
                                    A_domain_bottom << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                                    B_domain_bottom << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                                    C_domain_bottom << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                                    D_domain_bottom << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                                    itt[0] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom, tmp_A1);
                                    itt[1] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom, tmp_B1);
                                    //std::cout << "BOTTOM: " << itt[0] <<", " << itt[1] << std::endl <<std::endl;;
                                    if (itt[0] == true && itt[1] == true)
                                    {
                                        break;
                                    };
                                }

                                if (BC_DOMAIN_FLAG[2] == 1)
                                {
                                    Vector3d A_domain_front, B_domain_front, C_domain_front, D_domain_front;
                                    A_domain_front << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                                    B_domain_front << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                                    C_domain_front << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                                    D_domain_front << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                                    itt[0] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_front, B_domain_front, C_domain_front, D_domain_front, tmp_A1);
                                    itt[1] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_front, B_domain_front, C_domain_front, D_domain_front, tmp_B1);

                                    if (itt[0] == true && itt[1] == true)
                                    {
                                        break;
                                    };
                                }

                                if (BC_DOMAIN_FLAG[3] == 1)
                                {
                                    Vector3d A_domain_back, B_domain_back, C_domain_back, D_domain_back;
                                    A_domain_back << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                                    B_domain_back << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                                    C_domain_back << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                                    D_domain_back << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                                    itt[0] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_back, B_domain_back, C_domain_back, D_domain_back, tmp_A1);
                                    itt[1] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_back, B_domain_back, C_domain_back, D_domain_back, tmp_B1);

                                    if (itt[0] == true && itt[1] == true)
                                    {
                                        break;
                                    };
                                }

                                if (BC_DOMAIN_FLAG[4] == 1)
                                {
                                    Vector3d A_domain_left, B_domain_left, C_domain_left, D_domain_left;
                                    A_domain_left << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
                                    B_domain_left << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);
                                    C_domain_left << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
                                    D_domain_left << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
                                    itt[0] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_left, B_domain_left, C_domain_left, D_domain_left, tmp_A1);
                                    itt[1] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_left, B_domain_left, C_domain_left, D_domain_left, tmp_B1);

                                    if (itt[0] == true && itt[1] == true)
                                    {
                                        break;
                                    };
                                }

                                if (BC_DOMAIN_FLAG[5] == 1)
                                {
                                    Vector3d A_domain_right, B_domain_right, C_domain_right, D_domain_right;
                                    A_domain_right << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
                                    B_domain_right << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
                                    C_domain_right << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
                                    D_domain_right << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
                                    itt[0] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_right, B_domain_right, C_domain_right, D_domain_right, tmp_A1);
                                    itt[1] = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_right, B_domain_right, C_domain_right, D_domain_right, tmp_B1);

                                    if (itt[0] == true && itt[1] == true)
                                    {
                                        break;
                                    };
                                }
                            }

                            //-----------------------------------------------------------------------------------------------------
                            if (itt[0] == false && itt[1] == false)
                            {
                                Vector5d tmp_yu;
                                tmp_yu << 0, 0, 0, 0, 0;
                                JBA_2.insert(std::pair<std::pair<size_t, size_t>, Vector5d>(tmp_pair, tmp_yu));
                                //std::pair<std::map<std::pair<size_t, size_t>, Vector5d>::iterator, bool> ret = JBA_2.insert(std::pair<std::pair<size_t, size_t>, Vector5d>(tmp_pair, tmp_yu));
                                /*
                                    if (ret.second == 1)
                                    {
                                        std::cout << "the " << k << "th face, edge->" << EDGE_NO << ", element NO->" << j<< " : ";
                                        while (it_A != itEnd_A)
                                        {
                                            std::cout << it_A->first << ", [";
                                            std::cout << JXY_3D_1[JM_1[j](it_A->first)] << "], ";
                                            it_A++;
                                        }
                                        std::cout << "\n\n\n";
                                    }
                                    */
                            }
                        }
                    };
                };

                //fourth, traces' node
                ///< std::vector<std::pair<Vector3d, Vector3d>> Connection_traces;
                ///< std::vector<std::vector<size_t>> Trace_Node_A;
                Trace_Node_A.resize(Connection_traces.size());
                for (size_t j = 0; j < JXY_3D_1.size(); ++j)
                {
                    Vector3d PNT2D = JXY_1[j];
                    for (size_t k = 0; k < Connection_traces.size(); ++k)
                    {
                        bool qu = Determine_IF_a_pnt_lies_in_a_2D_line_segment(Connection_traces[k].first, Connection_traces[k].second, PNT2D);
                        if (qu == true)
                        {
                            Trace_Node_A[k].push_back(j);
                        }
                    }
                }

                //std::vector<Vector2d> JBA_1; ///< 1st kind of BC
                //std::vector<Vector6d> JBA_2; ///< 2nd kind of BC
                JXY.push_back(JXY_1);
                JXY_3D.push_back(JXY_3D_1);
                JM.push_back(JM_1);
                JB_1.push_back(JBA_1);
                JB_2.push_back(JBA_2);
                Trace_Node_overall.push_back(Trace_Node_A);
                Rota_angle.push_back(Rota_angle_1);
                neigh_shared.push_back(neigh_shared_A);
                NO_Nodes_p_each_frac.push_back(NO_Nodes_p_A);
            }
        }
    }

    Matlab_Plot_DFN_mesh("tdfn01_DFN_mesh.m", dom, JXY_3D, JM, JB_1, JB_2, Trace_Node_overall, P_idx);

    //void Address_Repetitive_Node(std::vector<std::vector<Vector3d>> &Coe_Matr_guide, const std::vector<std::vector<Vector3d>> JXY, size_t &overall_matrix_dimension);
    Address_Repetitive_Node(Coe_Matr_guide, Coe_Matr_guide_m, JXY_3D, NO_all_pnts);
    /*
        std::cout << "guide\n";
        for(size_t i = 0; i < Coe_Matr_guide.size(); ++i)
        {
            for(size_t j = 0; j < Coe_Matr_guide[i].size(); ++j)
            {
                std::cout << Coe_Matr_guide[i][j] << ", ";
            }
            std::cout << std::endl;
        }*/

    // now, identify boundary elements again, aiming at calculate inlet and outlei flux rates
    //std::vector<std::map<std::pair<size_t, size_t>, double>> Inlet; // the ith frac, and the jth element, the length of the edge
    //std::vector<std::map<std::pair<size_t, size_t>, double>> Outlet;

    for (size_t i = 0; i < JXY_3D.size(); ++i)
    {
        if (dom.Fractures[JXY_frac_ID[i]].If_intersect_surfaces(0) == 1) // top
        {
            std::map<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>> Inlet_UY;
            Vector3d A_domain_top, B_domain_top, C_domain_top, D_domain_top;
            A_domain_top << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(0);
            B_domain_top << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(0);
            C_domain_top << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(0);
            D_domain_top << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(0);
            for (size_t k = 0; k < JM[i].size(); ++k)
            {
                size_t IDX_count = 0;
                std::vector<size_t> ss;

                for (size_t j = 0; j < 6; j += 2) // 0 2 4
                {
                    Vector3d tmp_A1;
                    /*std::cout << "1.05\n";
                    std::cout << "frac id = " << i << ", size = " << JXY_3D[i].size() << std::endl;
                    std::cout << "frac id = " << i << ", node id = " << JM[i][k](j) << std::endl;
                    std::cout << "JM[i][k]: " << JM[i][k](0) << ", " << JM[i][k](2) << ", " << JM[i][k](4) << std::endl;*/
                    tmp_A1 = JXY_3D[i][JM[i][k](j)];
                    //std::cout << "1.1\n";
                    bool itt = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_top, B_domain_top, C_domain_top, D_domain_top, tmp_A1);
                    if (itt == true)
                    {
                        ss.push_back(j);
                        IDX_count++;
                    }
                }

                if (IDX_count == 2)
                {
                    Vector6d SEGMENT_2D;
                    Vector6d SEGMENT_3D;
                    //double edge_Length = 0;
                    if (ss[0] == 0 && ss[1] == 2)
                    {
                        SEGMENT_3D(0) = JXY_3D[i][JM[i][k](0)](0);
                        SEGMENT_3D(1) = JXY_3D[i][JM[i][k](0)](1);
                        SEGMENT_3D(2) = JXY_3D[i][JM[i][k](0)](2);
                        SEGMENT_3D(3) = JXY_3D[i][JM[i][k](2)](0);
                        SEGMENT_3D(4) = JXY_3D[i][JM[i][k](2)](1);
                        SEGMENT_3D(5) = JXY_3D[i][JM[i][k](2)](2);

                        SEGMENT_2D(0) = JXY[i][JM[i][k](0)](0);
                        SEGMENT_2D(1) = JXY[i][JM[i][k](0)](1);
                        SEGMENT_2D(2) = 0;
                        SEGMENT_2D(3) = JXY[i][JM[i][k](2)](0);
                        SEGMENT_2D(4) = JXY[i][JM[i][k](2)](1);
                        SEGMENT_2D(5) = JXY[i][JM[i][k](2)](2);

                        /*
                        Vector3d pnt3 = pnt2 - pnt1;
                        edge_Length = pow(pow(pnt3(0), 2) + pow(pnt3(1), 2) + pow(pnt3(2), 2), 0.5);
                        */
                    }
                    else if (ss[0] == 2 && ss[1] == 4)
                    {
                        /*
                        Vector3d pnt1 = JXY_3D[i][JM[i][k](4)];
                        Vector3d pnt2 = JXY_3D[i][JM[i][k](2)];*/

                        SEGMENT_3D(0) = JXY_3D[i][JM[i][k](4)](0);
                        SEGMENT_3D(1) = JXY_3D[i][JM[i][k](4)](1);
                        SEGMENT_3D(2) = JXY_3D[i][JM[i][k](4)](2);
                        SEGMENT_3D(3) = JXY_3D[i][JM[i][k](2)](0);
                        SEGMENT_3D(4) = JXY_3D[i][JM[i][k](2)](1);
                        SEGMENT_3D(5) = JXY_3D[i][JM[i][k](2)](2);

                        SEGMENT_2D(0) = JXY[i][JM[i][k](4)](0);
                        SEGMENT_2D(1) = JXY[i][JM[i][k](4)](1);
                        SEGMENT_2D(2) = 1;
                        SEGMENT_2D(3) = JXY[i][JM[i][k](2)](0);
                        SEGMENT_2D(4) = JXY[i][JM[i][k](2)](1);
                        SEGMENT_2D(5) = JXY[i][JM[i][k](2)](2);
                        /*
                        Vector3d pnt3 = pnt2 - pnt1;
                        edge_Length = pow(pow(pnt3(0), 2) + pow(pnt3(1), 2) + pow(pnt3(2), 2), 0.5);
                        */
                    }
                    else if (ss[0] == 0 && ss[1] == 4)
                    {
                        /*
                        Vector3d pnt1 = JXY_3D[i][JM[i][k](0)];
                        Vector3d pnt2 = JXY_3D[i][JM[i][k](4)];
                        Vector3d pnt3 = pnt2 - pnt1;
                        edge_Length = pow(pow(pnt3(0), 2) + pow(pnt3(1), 2) + pow(pnt3(2), 2), 0.5);*/
                        SEGMENT_3D(0) = JXY_3D[i][JM[i][k](4)](0);
                        SEGMENT_3D(1) = JXY_3D[i][JM[i][k](4)](1);
                        SEGMENT_3D(2) = JXY_3D[i][JM[i][k](4)](2);
                        SEGMENT_3D(3) = JXY_3D[i][JM[i][k](0)](0);
                        SEGMENT_3D(4) = JXY_3D[i][JM[i][k](0)](1);
                        SEGMENT_3D(5) = JXY_3D[i][JM[i][k](0)](2);

                        SEGMENT_2D(0) = JXY[i][JM[i][k](4)](0);
                        SEGMENT_2D(1) = JXY[i][JM[i][k](4)](1);
                        SEGMENT_2D(2) = 2;
                        SEGMENT_2D(3) = JXY[i][JM[i][k](0)](0);
                        SEGMENT_2D(4) = JXY[i][JM[i][k](0)](1);
                        SEGMENT_2D(5) = JXY[i][JM[i][k](0)](2);
                    }
                    //std::cout << "found a element top; element no: " << k << ";\n";
                    std::pair<size_t, size_t> tmp_pair1 = std::make_pair(i, k); // the ith frac, the k th element
                    std::pair<Vector6d, Vector6d> tmp_pair2 = std::make_pair(SEGMENT_2D, SEGMENT_3D);
                    Inlet_UY.insert(std::pair<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>>(tmp_pair1, tmp_pair2));
                }
                else if (IDX_count > 2)
                {
                    std::cout << "Error! It is impossible that more than two vertices of a triangular element are lying on an edge or a plane!\n";
                    exit(0);
                }
            }
            Inlet.push_back(Inlet_UY);
        }
        //std::cout << "111;\n";

        if (dom.Fractures[JXY_frac_ID[i]].If_intersect_surfaces(1) == 1) // bottom
        {
            std::map<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>> Outlet_UY;
            Vector3d A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom;
            A_domain_bottom << dom.Model_domain(4), dom.Model_domain(2), dom.Model_domain(1);
            B_domain_bottom << dom.Model_domain(5), dom.Model_domain(2), dom.Model_domain(1);
            C_domain_bottom << dom.Model_domain(5), dom.Model_domain(3), dom.Model_domain(1);
            D_domain_bottom << dom.Model_domain(4), dom.Model_domain(3), dom.Model_domain(1);

            for (size_t k = 0; k < JM[i].size(); ++k)
            {
                size_t IDX_count = 0;
                std::vector<size_t> ss;
                for (size_t j = 0; j < 6; j += 2)
                {

                    Vector3d tmp_A1;
                    tmp_A1 = JXY_3D[i][JM[i][k](j)];
                    bool itt = Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(A_domain_bottom, B_domain_bottom, C_domain_bottom, D_domain_bottom, tmp_A1);
                    if (itt == true)
                    {
                        ss.push_back(j);
                        IDX_count++;
                    }
                }

                if (IDX_count == 2)
                {
                    //std::cout << "found a element bottom; element no: " << k << ";\n";
                    //double edge_Length = 0;
                    Vector6d SEGMENT_2D;
                    Vector6d SEGMENT_3D;
                    if (ss[0] == 0 && ss[1] == 2)
                    {
                        SEGMENT_3D(0) = JXY_3D[i][JM[i][k](0)](0);
                        SEGMENT_3D(1) = JXY_3D[i][JM[i][k](0)](1);
                        SEGMENT_3D(2) = JXY_3D[i][JM[i][k](0)](2);
                        SEGMENT_3D(3) = JXY_3D[i][JM[i][k](2)](0);
                        SEGMENT_3D(4) = JXY_3D[i][JM[i][k](2)](1);
                        SEGMENT_3D(5) = JXY_3D[i][JM[i][k](2)](2);

                        SEGMENT_2D(0) = JXY[i][JM[i][k](0)](0);
                        SEGMENT_2D(1) = JXY[i][JM[i][k](0)](1);
                        SEGMENT_2D(2) = 0;
                        SEGMENT_2D(3) = JXY[i][JM[i][k](2)](0);
                        SEGMENT_2D(4) = JXY[i][JM[i][k](2)](1);
                        SEGMENT_2D(5) = JXY[i][JM[i][k](2)](2);
                        /*
                        Vector3d pnt1 = JXY_3D[i][JM[i][k](0)];
                        Vector3d pnt2 = JXY_3D[i][JM[i][k](2)];
                        Vector3d pnt3 = pnt2 - pnt1;
                        edge_Length = pow(pow(pnt3(0), 2) + pow(pnt3(1), 2) + pow(pnt3(2), 2), 0.5);*/
                    }
                    else if (ss[0] == 2 && ss[1] == 4)
                    {
                        SEGMENT_3D(0) = JXY_3D[i][JM[i][k](4)](0);
                        SEGMENT_3D(1) = JXY_3D[i][JM[i][k](4)](1);
                        SEGMENT_3D(2) = JXY_3D[i][JM[i][k](4)](2);
                        SEGMENT_3D(3) = JXY_3D[i][JM[i][k](2)](0);
                        SEGMENT_3D(4) = JXY_3D[i][JM[i][k](2)](1);
                        SEGMENT_3D(5) = JXY_3D[i][JM[i][k](2)](2);

                        SEGMENT_2D(0) = JXY[i][JM[i][k](4)](0);
                        SEGMENT_2D(1) = JXY[i][JM[i][k](4)](1);
                        SEGMENT_2D(2) = 1;
                        SEGMENT_2D(3) = JXY[i][JM[i][k](2)](0);
                        SEGMENT_2D(4) = JXY[i][JM[i][k](2)](1);
                        SEGMENT_2D(5) = JXY[i][JM[i][k](2)](2);
                        /*
                        Vector3d pnt1 = JXY_3D[i][JM[i][k](4)];
                        Vector3d pnt2 = JXY_3D[i][JM[i][k](2)];
                        Vector3d pnt3 = pnt2 - pnt1;
                        edge_Length = pow(pow(pnt3(0), 2) + pow(pnt3(1), 2) + pow(pnt3(2), 2), 0.5);*/
                    }
                    else if (ss[0] == 0 && ss[1] == 4)
                    {
                        SEGMENT_3D(0) = JXY_3D[i][JM[i][k](4)](0);
                        SEGMENT_3D(1) = JXY_3D[i][JM[i][k](4)](1);
                        SEGMENT_3D(2) = JXY_3D[i][JM[i][k](4)](2);
                        SEGMENT_3D(3) = JXY_3D[i][JM[i][k](0)](0);
                        SEGMENT_3D(4) = JXY_3D[i][JM[i][k](0)](1);
                        SEGMENT_3D(5) = JXY_3D[i][JM[i][k](0)](2);

                        SEGMENT_2D(0) = JXY[i][JM[i][k](4)](0);
                        SEGMENT_2D(1) = JXY[i][JM[i][k](4)](1);
                        SEGMENT_2D(2) = 2;
                        SEGMENT_2D(3) = JXY[i][JM[i][k](0)](0);
                        SEGMENT_2D(4) = JXY[i][JM[i][k](0)](1);
                        SEGMENT_2D(5) = JXY[i][JM[i][k](0)](2);
                        /*
                        Vector3d pnt1 = JXY_3D[i][JM[i][k](0)];
                        Vector3d pnt2 = JXY_3D[i][JM[i][k](4)];
                        Vector3d pnt3 = pnt2 - pnt1;
                        edge_Length = pow(pow(pnt3(0), 2) + pow(pnt3(1), 2) + pow(pnt3(2), 2), 0.5);*/
                    }
                    std::pair<size_t, size_t> tmp_pair1 = std::make_pair(i, k); // the ith frac, the k th element
                    std::pair<Vector6d, Vector6d> tmp_pair2 = std::make_pair(SEGMENT_2D, SEGMENT_3D);
                    Outlet_UY.insert(std::pair<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>>(tmp_pair1, tmp_pair2));
                }
                else if (IDX_count > 2)
                {
                    std::cout << "Error! It is impossible that more than two vertices of a triangular element are lying on an edge or a plane!\n";
                    exit(0);
                }
            }
            Outlet.push_back(Outlet_UY);
        }
        //std::cout << "222;\n";
    }
    //std::cout << "333;\n";

    /*
        std::cout << "top: " << Inlet.size() << std::endl;
        std::map<std::pair<size_t, size_t>, double>::iterator iter;
        iter = Inlet.begin();
        while (iter != Inlet.end())
        {
            std::cout << iter->first.first << ", " << iter->first.second << ";\n";
            iter++;
        }
        std::cout << "bottom: " << Outlet.size() << std::endl;*/
};
///************************************************************///

inline void DFN_mesh::Find_pnt_inside_region(const DFN::Domain dom, const std::vector<Vector3d> tem_verts_trim, const size_t i, const size_t P_idx, const size_t he, double &x_ik, double &y_ik)
{
    Vector3d Center_p;
    Center_p << 0, 0, 0;
    if (pointInConvexPolygon(Center_p, tem_verts_trim) == true)
    {
        x_ik = 0;
        y_ik = 0;
        return;
    }

    if (dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Verts_trim.size() == 3)
    {
        Vector3d A, B, C, D;
        A << tem_verts_trim[0](0), tem_verts_trim[0](1), tem_verts_trim[0](2);
        B << (tem_verts_trim[1](0) + tem_verts_trim[2](0)) / 0.5, (tem_verts_trim[1](1) + tem_verts_trim[2](1)) / 0.5, (tem_verts_trim[1](2) + tem_verts_trim[2](2)) / 0.5;
        C << tem_verts_trim[1](0), tem_verts_trim[1](1), tem_verts_trim[1](2);
        D << (tem_verts_trim[0](0) + tem_verts_trim[2](0)) / 0.5, (tem_verts_trim[0](1) + tem_verts_trim[2](1)) / 0.5, (tem_verts_trim[0](2) + tem_verts_trim[2](2)) / 0.5;

        Line A1(A, B);
        Line A2(C, D);
        size_t yi = is_intersect(A1, A2);
        if (yi == 1)
            intersection_between_two_line_segments(A(0), A(1), B(0), B(1), C(0), C(1), D(0), D(1), x_ik, y_ik);
        else
        {
            std::cout << "Mesh_Frac--> cannot find a point inside the fracture\n";
            exit(0);
        }
    }
    else
    {
        size_t j = 1;
    gkp0:;
        if (j > dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Verts_trim.size() - 1 - 2)
        {
            std::cout << "Mesh_Frac--> cannot find a point inside the fracture\n";
            exit(0);
        }
        Vector3d A, B, C, D;
        A << tem_verts_trim[0](0), tem_verts_trim[0](1), tem_verts_trim[0](2);
        B << tem_verts_trim[2](0), tem_verts_trim[2](1), tem_verts_trim[2](2);
        C << tem_verts_trim[j](0), tem_verts_trim[j](1), tem_verts_trim[j](2);
        D << tem_verts_trim[j + 2](0), tem_verts_trim[j + 2](1), tem_verts_trim[j + 2](2);

        Line A1(A, B);
        Line A2(C, D);
        size_t yi = is_intersect(A1, A2);
        if (yi == 1)
            intersection_between_two_line_segments(A(0), A(1), B(0), B(1), C(0), C(1), D(0), D(1), x_ik, y_ik);
        else
        {
            j++;
            goto gkp0;
        }
    }
};

inline void DFN_mesh::Rotation_3D_frac_to_2D(DFN::Domain dom, std::vector<Vector3d> &tem_verts_trim_AAA, std::map<std::pair<double, double>, int> &MapPnt, std::vector<Vector2d> &SegIDtoTD, Vector3d &temp3, Vector3d &Normal_frac, const size_t i, const size_t P_idx, const size_t he, const double d_In, std::vector<std::pair<Vector3d, Vector3d>> &Connection_traces)
{
    std::vector<Vector3d> tem_verts_trim;
    std::vector<std::vector<Vector3d>> Seg_pnt; ///< segment ends (not include domain coordinates)
    std::vector<Vector3d> Iso_pnts;
    size_t NoIT = 0; ///< number_of_intersection_traces
    tem_verts_trim.resize(dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Verts_trim.size());
    std::vector<int> no_overlaped_pnt;
    Normal_frac = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Normal_vector;
    Find_vector_2(dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Normal_vector, temp3);
    for (size_t j = 0; j < dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Verts_trim.size(); ++j)
    {
        tem_verts_trim[j] = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Verts_trim[j] - dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center;
        if ((j > 0) && (abs(tem_verts_trim[j](0) - tem_verts_trim[j - 1](0)) < 0.0001) && (abs(tem_verts_trim[j](1) - tem_verts_trim[j - 1](1)) < 0.0001))
        {
            //nothing
        }
        else
        {
            no_overlaped_pnt.push_back(j);
        }
    }
    tem_verts_trim_AAA.resize(no_overlaped_pnt.size());
    //we first need to know if all the traces intersect or not
    //for each trace, how many truncation points there are
    std::vector<std::vector<Vector3d>> Intersec_traces; ///< the first index is trace ID, and the second is the intersections between this trace and other traces (if they intersects)
    std::vector<std::pair<Vector3d, Vector3d>> Traces;  ///< the two coordinates of traces
    //so, a intersection trace contains points: Traces.first, Intersec_traces (point), Traces.second
    //we need a Tag to indicate the attribute of each point; -2 means they are edge ends, but they are really the same (i.e., intersection is a point)
    //-1 means edge ends, not the same
    // 0 just means the common intersection points between traces
    std::vector<std::vector<int>> Tag_attri;
    std::vector<int> Tag_edge;
    //std::cout << "Trace ends: \n";

    for (size_t j = 0; j < dom.Connections.size(); j += 2)
    {
        if (dom.Connections[j] == dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Tag || dom.Connections[j + 1] == dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Tag)
        {
            NoIT++;
            std::pair<size_t, size_t> sp = std::make_pair(dom.Connections[j], dom.Connections[j + 1]);
            // dom.Intersections[sp].first, dom.Intersections[sp].second
            if ((abs(dom.Intersections[sp].first(0) - dom.Intersections[sp].second(0)) < 0.0001) && (abs(dom.Intersections[sp].first(1) - dom.Intersections[sp].second(1)) < 0.0001))
            {
                std::pair<Vector3d, Vector3d> trace_tmp = std::make_pair(dom.Intersections[sp].first, dom.Intersections[sp].second);
                Traces.push_back(trace_tmp);
                Tag_edge.push_back(1); // means this edge is a point
                //std::cout << dom.Intersections[sp].first << ", " << dom.Intersections[sp].second << "\n";
            }
            else
            {
                std::pair<Vector3d, Vector3d> trace_tmp = std::make_pair(dom.Intersections[sp].first, dom.Intersections[sp].second);
                Traces.push_back(trace_tmp);
                Tag_edge.push_back(0);
                //std::cout << dom.Intersections[sp].first << ", " << dom.Intersections[sp].second << "\n";
            }
            std::pair<Vector3d, Vector3d> Ap = std::make_pair(dom.Intersections[sp].first, dom.Intersections[sp].second);
            Connection_traces.push_back(Ap);
        }
    }

    //from 3d to 2d
    if (abs(Normal_frac(0)) < 0.000001 && abs(Normal_frac(1)) < 0.000001 && abs(Normal_frac(2)) < 0.000001)
    {
        //Verts;
    }
    else
    {
        double R_angle_temp1 = 0;
        double x_temp = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Dip_angle; ///it is better to create a new variable to represent the dip angle, because debuging shows direct use of 'Dip_angle' to calculate rotation angle leads wrong output
        R_angle_temp1 = -x_temp * M_PI / 180;

        Quaternion_t Q_axis_1;
        NormalizeRotation(R_angle_temp1, temp3, Q_axis_1);

        for (size_t j = 0; j < Traces.size(); j++)
        {
            Vector3d temp4;
            Rotation(Traces[j].first - dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center, Q_axis_1, temp4);
            Traces[j].first = temp4;

            Vector3d temp5;
            Rotation(Traces[j].second - dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center, Q_axis_1, temp5);
            Traces[j].second = temp5;
        };

        for (size_t j = 0; j < tem_verts_trim.size(); j++)
        {
            Vector3d temp4;
            Rotation(tem_verts_trim[j], Q_axis_1, temp4);
            tem_verts_trim[j] = temp4;
            //tem_verts_trim_AAA[j] = tem_verts_trim[j];
        };

        for (size_t j = 0; j < Connection_traces.size(); ++j)
        {
            Vector3d temp4;
            Rotation(Connection_traces[j].first - dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center, Q_axis_1, temp4);
            Connection_traces[j].first = temp4;

            Vector3d temp5;
            Rotation(Connection_traces[j].second - dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][i]].Center, Q_axis_1, temp5);
            Connection_traces[j].second = temp5;
        }
    }

    for (size_t j = 0; j < tem_verts_trim_AAA.size(); ++j)
    {
        tem_verts_trim_AAA[j] = tem_verts_trim[no_overlaped_pnt[j]];
    }

    //std::vector<std::pair<Vector3d, Vector3d>> Traces; ///< the two coordinates of traces
    std::vector<Vector3d> tem_verts_trim_BBB;
    for (size_t j = 0; j < tem_verts_trim_AAA.size(); ++j)
    {
        tem_verts_trim_BBB.push_back(tem_verts_trim_AAA[j]);

        int j_plus_1 = (j + 1) % tem_verts_trim_AAA.size();
        //if bound intersects traces
        std::vector<Vector3d> tmp_intersect_bound_trace;
        //std::cout << "bound ID: " << j << std::endl;
        for (size_t k = 0; k < Traces.size(); ++k)
        {
            //std::cout << "Traces ID: " << k << std::endl;
            //std::cout << "bound ID: " << j << ", " <<j_plus_1  << std::endl;
            Line A_tmp_1(Traces[k].first, Traces[k].second);
            Line B_tmp_1(tem_verts_trim_AAA[j], tem_verts_trim_AAA[j_plus_1]);
            size_t tmp_yi = is_intersect(A_tmp_1, B_tmp_1);
            //std::cout << "tmp_yi: " << tmp_yi <<std::endl;
            if (tmp_yi == 1)
            {
                //std::cout<<"find intersection between bound and trace!\n";
                double x_yu_1, y_yu_1;
                intersection_between_two_line_segments(A_tmp_1.xa, A_tmp_1.ya, A_tmp_1.xb, A_tmp_1.yb, B_tmp_1.xa, B_tmp_1.ya, B_tmp_1.xb, B_tmp_1.yb, x_yu_1, y_yu_1);
                Vector3d tmp_ty_2;
                tmp_ty_2 << x_yu_1, y_yu_1, 0;
                tmp_intersect_bound_trace.push_back(tmp_ty_2);
                //std::cout << tmp_ty_2 << std::endl;
            }
        }
        Sort_trace_intersection_pnts(tmp_intersect_bound_trace, tem_verts_trim_AAA[j]);
        //std::cout << "tmp_intersect_bound_trace.size(): " << tmp_intersect_bound_trace.size() << std::endl;
        for (size_t k = 0; k < tmp_intersect_bound_trace.size(); ++k)
        {
            tem_verts_trim_BBB.push_back(tmp_intersect_bound_trace[k]);
        }
    }
    //std::cout << "tem_verts_trim_BBB.size(): " << tem_verts_trim_BBB.size() << std::endl;

    //we need to identify if there are intersections between domain bounds and fracture intersection traces
    tem_verts_trim.resize(0);
    for (size_t j = 0; j < tem_verts_trim_BBB.size(); j++)
    {
        Vector3d initial_pnt;
        Vector3d end_pnt;

        initial_pnt = tem_verts_trim_BBB[j];
        end_pnt = tem_verts_trim_BBB[(j + 1) % tem_verts_trim_BBB.size()];
        tem_verts_trim.push_back(tem_verts_trim_BBB[j]); // push initial point

        double len_AAA = Length_2d_segment(initial_pnt, end_pnt);
        int N_AAA = len_AAA / d_In + 1.;
        std::vector<Vector3d> tmp_AAA;
        binary_linear_equation(initial_pnt, end_pnt, N_AAA, tmp_AAA);

        for (size_t k = 0; k < tmp_AAA.size(); ++k)
            tem_verts_trim.push_back(tmp_AAA[k]);
    }

    //std::cout << "NoIT: " << NoIT << std::endl;
    Intersec_traces.resize(NoIT);
    for (size_t j = 0; j < NoIT; ++j)
    {
        for (size_t k = 0; k < NoIT; ++k)
        {
            if (j != k)
            {
                Line A_tmp_1(Traces[j].first, Traces[j].second);
                Line B_tmp_1(Traces[k].first, Traces[k].second);
                size_t tmp_yi = is_intersect(A_tmp_1, B_tmp_1);

                if (tmp_yi == 1)
                {
                    if ((abs(Traces[j].first(0) - Traces[j].second(0)) < 0.0001 && abs(Traces[j].first(1) - Traces[j].second(1)) < 0.0001) && (abs(Traces[k].first(0) - Traces[k].second(0)) > 0.0001 || abs(Traces[k].first(1) - Traces[k].second(1)) > 0.0001))
                    {
                        //if trace J is a point
                        Intersec_traces[j].push_back(Traces[j].first);
                    }
                    else if ((abs(Traces[k].first(0) - Traces[k].second(0)) < 0.0001 && abs(Traces[k].first(1) - Traces[k].second(1)) < 0.0001) && (abs(Traces[j].first(0) - Traces[j].second(0)) > 0.0001 || abs(Traces[j].first(1) - Traces[j].second(1)) > 0.0001))
                    {
                        //if trace k is a point
                        Intersec_traces[j].push_back(Traces[k].first);
                    }
                    else if ((abs(Traces[j].first(0) - Traces[j].second(0)) < 0.0001 && abs(Traces[j].first(1) - Traces[j].second(1)) < 0.0001) && (abs(Traces[k].first(0) - Traces[k].second(0)) < 0.0001 && abs(Traces[k].first(1) - Traces[k].second(1)) < 0.0001))
                    {
                        //if J and K are both points
                        Intersec_traces[j].push_back(Traces[k].first);
                    }
                    else
                    {
                        //std::cout << "intersection between traces found!\n";
                        double x_yu_1, y_yu_1;
                        intersection_between_two_line_segments(Traces[j].first(0), Traces[j].first(1), Traces[j].second(0), Traces[j].second(1), Traces[k].first(0), Traces[k].first(1), Traces[k].second(0), Traces[k].second(1), x_yu_1, y_yu_1);
                        Vector3d tmp_ty_2;
                        tmp_ty_2 << x_yu_1, y_yu_1, 0;
                        Intersec_traces[j].push_back(tmp_ty_2);
                        //std::cout << tmp_ty_2 << std::endl;
                    }
                }
            }
        }
    }

    //sort
    for (size_t j = 0; j < Intersec_traces.size(); ++j)
    {
        Sort_trace_intersection_pnts(Intersec_traces[j], Traces[j].first);
    }

    Tag_attri.resize(Intersec_traces.size());
    for (size_t j = 0; j < Tag_attri.size(); ++j)
    {
        Tag_attri[j].resize(2. + Intersec_traces[j].size());
        if (Tag_edge[j] == 0) // this trace is a segment
        {
            Tag_attri[j][0] = -1;                       //edge ends, not the same
            Tag_attri[j][Tag_attri[j].size() - 1] = -1; //edge ends, not the same
        }
        else //this trace is a point
        {
            Tag_attri[j][0] = -2; // -2 means a point
            Tag_attri[j][Tag_attri[j].size() - 1] = -2;
        }
    }

    /// and then push them into matrix of vertices
    // ----- std::vector<std::vector<Vector3d>> Intersec_traces;    ///< the first index is trace ID, and the second is the intersections between this trace and other traces (if they intersects)
    // ----- std::vector<std::pair<Vector3d, Vector3d>> Traces; ///< the two coordinates of traces
    // R0tation Intersec_traces and Traces
    Seg_pnt.resize(NoIT);
    for (size_t j = 0; j < NoIT; ++j) // NoIT: number of traces
    {
        //tem_verts_trim.push_back(dom.Intersections[sp].first - dom.Fractures[i].Center)
        if (Tag_attri[j][0] == -1) // it is a segment but not a domain boundary
        {
            Seg_pnt[j].push_back(Traces[j].first);
            //std::cout << " A = [" << Traces[j].first(0) << ", "<< Traces[j].first(1) << "; \n";
            //Seg_pnt.push_back(Traces[j].second);
            if (Intersec_traces[j].size() != 0)
            {
                for (size_t k = 0; k < Intersec_traces[j].size(); ++k)
                {

                    if (k == 0)
                    {
                        Vector3d tmp_A;
                        std::vector<Vector3d> tmp_Inserted;
                        tmp_A = Intersec_traces[j][k] - Traces[j].first;
                        double length_a = pow(pow(tmp_A(0), 2.) + pow(tmp_A(1), 2.), .5);
                        int a = (int)(length_a / d_In);
                        a++; ///< a is the number of subsegments along a intersection traces; when a  == 1, means no segment

                        binary_linear_equation(Traces[j].first, Intersec_traces[j][k], a, tmp_Inserted);
                        for (size_t ko = 0; ko < tmp_Inserted.size(); ++ko)
                        {
                            Seg_pnt[j].push_back(tmp_Inserted[ko]);
                            //std::cout << tmp_Inserted[ko](0) << ", " << tmp_Inserted[ko](1) << ";\n";
                        }
                        Seg_pnt[j].push_back(Intersec_traces[j][k]);
                        //std::cout << Intersec_traces[j][k](0) << ", " << Intersec_traces[j][k](1) << ";\n";
                        //std::cout << "splitting: " << Traces[j].first << ", " << Intersec_traces[j][k] << std::endl;
                    }
                    else if (k > 0 && k <= Intersec_traces[j].size() - 1)
                    {
                        Vector3d tmp_A;
                        std::vector<Vector3d> tmp_Inserted;
                        tmp_A = Intersec_traces[j][k] - Intersec_traces[j][k - 1];
                        double length_a = pow(pow(tmp_A(0), 2.) + pow(tmp_A(1), 2.), .5);
                        int a = (int)(length_a / d_In);
                        a++; ///< a is the number of subsegments along a intersection traces; when a  == 1, means no segment

                        binary_linear_equation((Intersec_traces[j][k - 1]), (Intersec_traces[j][k]), a, tmp_Inserted);
                        for (size_t ko = 0; ko < tmp_Inserted.size(); ++ko)
                        {
                            Seg_pnt[j].push_back(tmp_Inserted[ko]);
                            //std::cout << tmp_Inserted[ko](0) << ", " << tmp_Inserted[ko](1) << ";\n";
                        }
                        Seg_pnt[j].push_back(Intersec_traces[j][k]);
                        //std::cout << Intersec_traces[j][k](0) << ", " << Intersec_traces[j][k](1) << ";\n";
                        //std::cout << "splitting: " << (Intersec_traces[j][k - 1]) << ", " << (Intersec_traces[j][k]) << std::endl;
                    }
                }

                Vector3d tmp_A;
                std::vector<Vector3d> tmp_Inserted;
                tmp_A = Traces[j].second - Intersec_traces[j][Intersec_traces[j].size() - 1];
                double length_a = pow(pow(tmp_A(0), 2.) + pow(tmp_A(1), 2.), .5);
                int a = (int)(length_a / d_In);
                a++; ///< a is the number of subsegments along a intersection traces; when a  == 1, means no segment

                binary_linear_equation((Intersec_traces[j][Intersec_traces[j].size() - 1]), (Traces[j].second), a, tmp_Inserted);
                for (size_t ko = 0; ko < tmp_Inserted.size(); ++ko)
                {
                    Seg_pnt[j].push_back(tmp_Inserted[ko]);
                    //std::cout << tmp_Inserted[ko](0) << ", " << tmp_Inserted[ko](1) << ";\n";
                }
                //std::cout << "splitting: " << (Intersec_traces[j][Intersec_traces[j].size() - 1]) << ", " << (Traces[j].second) << std::endl;
            }
            else
            {
                Vector3d tmp_A;
                std::vector<Vector3d> tmp_Inserted;
                tmp_A = Traces[j].second - Traces[j].first;
                double length_a = pow(pow(tmp_A(0), 2.) + pow(tmp_A(1), 2.), .5);
                //std::cout << " length_a = " << length_a << std::endl;
                int a = (int)(length_a / d_In);
                a++; ///< a is the number of subsegments along a intersection traces; when a  == 1, means no segment

                //std::cout << "a(number of subsegments)is " << a << std::endl;
                binary_linear_equation((Traces[j].first), (Traces[j].second), a, tmp_Inserted);
                for (size_t ko = 0; ko < tmp_Inserted.size(); ++ko)
                {
                    Seg_pnt[j].push_back(tmp_Inserted[ko]);
                }
            }
            Seg_pnt[j].push_back(Traces[j].second);
            //std::cout << Traces[j].second(0) << ", " << Traces[j].second(1) << "]\n";
            //exit(0);
        }
        else if (Tag_attri[j][0] == -2)
        {
            Iso_pnts.push_back(Traces[j].first);
        }
        else
        {
            std::cout << "Error NO. 2! Something weried happened! In Class - 'Mesh_frac', Function - 'Rotation_3D_frac_to_2D'!\n";
            std::cout << "Two edge ends have not been assigned attributes yet!\n";
            exit(0);
        }
    }

    int counter_pnt_ID = 0;
    for (size_t i = 0; i < tem_verts_trim.size(); ++i)
    {
        std::pair<double, double> tmp_pnt_1 = std::make_pair(round(tem_verts_trim[i](0), 4), round(tem_verts_trim[i](1), 4));
        MapPnt.insert(std::pair<std::pair<double, double>, int>(tmp_pnt_1, i));

        Vector2d tmp_IDtoID;
        tmp_IDtoID << (int)(i), (int)((i + 1) % tem_verts_trim.size());
        SegIDtoTD.push_back(tmp_IDtoID);
        if (i == tem_verts_trim.size() - 1)
        {
            counter_pnt_ID = i + 1;
        }
    }

    for (size_t i = 0; i < Seg_pnt.size(); ++i)
    {
        std::vector<int> tmp_segID;
        int tmp_j = 0;
        for (size_t ut = 0; ut < i; ut++)
        {
            tmp_j = tmp_j + Seg_pnt[ut].size();
        }

        for (size_t j = 0; j < Seg_pnt[i].size(); ++j)
        {
            std::pair<double, double> tmp_pnt_1 = std::make_pair(round(Seg_pnt[i][j](0), 4), round(Seg_pnt[i][j](1), 4));
            std::pair<std::map<std::pair<double, double>, int>::iterator, bool> ret = MapPnt.insert(std::pair<std::pair<double, double>, int>(tmp_pnt_1, tmp_j + j + counter_pnt_ID));

            if (!ret.second) // fail to insert
            {
                counter_pnt_ID--;
                std::map<std::pair<double, double>, int>::iterator iter;
                iter = MapPnt.find(tmp_pnt_1);
                if (iter != MapPnt.end())
                {
                    //the point ID is
                    tmp_segID.push_back(iter->second);
                }
                else
                {
                    std::cout << "Error happened in Class 'Mesh_frac', Function 'Rotation_3D_frac_to_2D'!\n";
                    exit(0);
                }
            }
            else //insert successfully
            {
                tmp_segID.push_back(tmp_j + j + counter_pnt_ID);
            }
            if (i == Seg_pnt.size() - 1 && j == Seg_pnt[i].size() - 1)
            {
                counter_pnt_ID = tmp_j + j + counter_pnt_ID + 1;
            }
        }
        for (size_t j = 0; j < tmp_segID.size() - 1; ++j)
        {
            Vector2d tmp_IDtoID;
            tmp_IDtoID << tmp_segID[j], tmp_segID[j + 1];
            SegIDtoTD.push_back(tmp_IDtoID);
        }
    }

    for (size_t i = 0; i < Iso_pnts.size(); ++i)
    {
        std::pair<double, double> tmp_pnt_1 = std::make_pair(round(Iso_pnts[i](0), 4), round(Iso_pnts[i](1), 4));
        std::pair<std::map<std::pair<double, double>, int>::iterator, bool> ret = MapPnt.insert(std::pair<std::pair<double, double>, int>(tmp_pnt_1, i + counter_pnt_ID));
        if (!ret.second) // fail to insert
        {
            counter_pnt_ID--;
        }
        else //insert successfully
        {
            //do nothing
        }
    }
    //Matlab_plot_2D_frac("tdfn_PNT", i, tem_verts_trim, Seg_pnt);
};

//below output matlab script to graph mesh results
void DFN_mesh::Matlab_Plot_DFN_mesh(string FileKey, const DFN::Domain dom, const std::vector<std::vector<Vector3d>> JXY, const std::vector<std::vector<Vector6d>> JM, std::vector<std::map<size_t, double>> JB_1, std::vector<std::map<std::pair<size_t, size_t>, Vector5d>> JB_2, const std::vector<std::vector<std::vector<size_t>>> Trace_Node_overall, const size_t P_idx)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    //Plotting the fractures

    for (size_t he = 0; he < dom.Percolation_cluster[P_idx].size(); he++)
    {
        for (size_t nf = 0; nf < dom.Listofclusters[dom.Percolation_cluster[P_idx][he]].size(); nf++)
        {
            if (Listofclusters_mesh_only[dom.Percolation_cluster[P_idx][he]][nf] != -1)
            {
                size_t n_verts = dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][nf]].Verts_trim.size();
                oss << "fill3([";
                for (size_t nv = 0; nv < n_verts + 1; ++nv)
                {
                    size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                    oss << dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][nf]].Verts_trim[nv_1](0) << " ";
                }
                oss << "],[";
                for (size_t nv = 0; nv < n_verts + 1; ++nv)
                {
                    size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                    oss << dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][nf]].Verts_trim[nv_1](1) << " ";
                }
                oss << "],[";
                for (size_t nv = 0; nv < n_verts + 1; ++nv)
                {
                    size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                    oss << dom.Fractures[dom.Listofclusters[dom.Percolation_cluster[P_idx][he]][nf]].Verts_trim[nv_1](2) << " ";
                }
                oss << "],[rand rand rand]);\ngrid on; hold on;\n";
            }
        }
    }

    //Plotting the model domain
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << dom.Surfaces[i].Verts[j](0) << " " << dom.Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](1) << " " << dom.Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](2) << " " << dom.Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on; hold on;\n";
        }
    }
    double xmin_1 = dom.Model_domain(4), xmax_1 = dom.Model_domain(5);
    double ymin_1 = dom.Model_domain(2), ymax_1 = dom.Model_domain(3);
    double zmin_1 = dom.Model_domain(1), zmax_1 = dom.Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

    //Plotting mesh (DFN)
    oss << "\n\n%***************mesh*****************\n";
    //std::cout << "jjj1;\n";
    for (size_t k = 0; k < JM.size(); ++k)
    {
        oss << "JM" << k + 1 << " = [";
        for (size_t i = 0; i < JM[k].size(); ++i)
        {
            for (size_t j = 0; j < 6; ++j)
            {
                oss << JM[k][i](j) + 1 << " ";
            }
            oss << ";";
        }
        oss << "];\n";
    }
    //std::cout << "jjj2;\n";
    for (size_t k = 0; k < JXY.size(); ++k)
    {
        oss << "JXY" << k + 1 << " = [";
        for (size_t i = 0; i < JXY[k].size(); ++i)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                oss << JXY[k][i](j) << " ";
            }
            oss << ";";
        }
        oss << "];\n";
    }
    //std::cout << "kkk1;\n";
    for (size_t k = 0; k < JM.size(); ++k)
    {
        oss << "[me" << k + 1 << ", ne" << k + 1 << "] = size(JM" << k + 1 << ");\n";
        oss << "AAA" << k + 1 << " = [rand rand rand]\n";
        oss << "for i = 1:me" << k + 1 << "\n";
        oss << "\tfor j = 1:ne" << k + 1 << "\n";
        oss << "\t\tx1=JXY" << k + 1 << "(JM" << k + 1 << "(i, j), 1)\n";
        oss << "\t\tx2=JXY" << k + 1 << "(JM" << k + 1 << "(i, (j + 1) - fix((j+1)/(ne" << k + 1 << "+1))*(ne" << k + 1 << ")), 1)\n";
        oss << "\t\ty1=JXY" << k + 1 << "(JM" << k + 1 << "(i, j), 2)\n";
        oss << "\t\ty2=JXY" << k + 1 << "(JM" << k + 1 << "(i, (j + 1) - fix((j+1)/(ne" << k + 1 << "+1))*(ne" << k + 1 << ")), 2)\n";
        oss << "\t\tz1=JXY" << k + 1 << "(JM" << k + 1 << "(i, j), 3)\n";
        oss << "\t\tz2=JXY" << k + 1 << "(JM" << k + 1 << "(i, (j + 1) - fix((j+1)/(ne" << k + 1 << "+1))*(ne" << k + 1 << ")), 3)\n";
        oss << "\t\tplot3([x1 x2],[y1 y2], [z1 z2],'color', AAA" << k + 1 << ",'Linewidth',0.05)\n";
        oss << "\t\thold on\n";
        oss << "\tend;\n";
        oss << "end;\n";
        //oss << "clear JXY" << k + 1 << " JM" << k + 1 << " m" << k + 1 << " n" << k + 1 << " AAA" << k + 1 << ";\n\n";
    }
    //std::cout << "kk2;\n";
    oss << "\n\n%***************boundary*****************\n";
    for (size_t k = 0; k < JB_1.size(); ++k)
    {
        oss << "JB_1st" << k + 1 << " = [";
        std::map<size_t, double>::iterator it_1;
        std::map<size_t, double>::iterator itEnd_1;
        it_1 = JB_1[k].begin();
        itEnd_1 = JB_1[k].end();
        while (it_1 != itEnd_1)
        {
            oss << it_1->first + 1 << ", " << it_1->second << "; ";
            it_1++;
        }
        oss << "];\n";

        oss << "[m" << k + 1 << ", n" << k + 1 << "] = size(JB_1st" << k + 1 << ");\n";
        oss << "for i = 1:m" << k + 1 << "\n";
        oss << "\txa1=JXY" << k + 1 << "(JB_1st" << k + 1 << "(i, 1), 1)\n";
        oss << "\tya1=JXY" << k + 1 << "(JB_1st" << k + 1 << "(i, 1), 2)\n";
        oss << "\tza1=JXY" << k + 1 << "(JB_1st" << k + 1 << "(i, 1), 3)\n";
        oss << "\tscatter3(xa1, ya1, za1, 30, 'k', 'filled');\n";
        oss << "\thold on;\n";
        oss << "end;\n";
    }
    //std::cout << "kkk3;\n";
    oss << "\nclear x1 x2 xa1 y1 y2 ya1 z1 z2 za1;\n\n";
    //std::vector<std::map<std::pair<size_t, size_t>, Vector4d>> JB_2;
    for (size_t k = 0; k < JB_2.size(); ++k)
    {
        oss << "JB_2st" << k + 1 << " = [\n";

        std::map<std::pair<size_t, size_t>, Vector5d>::iterator it_1;
        std::map<std::pair<size_t, size_t>, Vector5d>::iterator itEnd_1;
        it_1 = JB_2[k].begin();
        itEnd_1 = JB_2[k].end();
        while (it_1 != itEnd_1)
        {
            size_t tmp_ui = 0;
            if (it_1->first.second == 0)
            {
                tmp_ui = 0;
            }
            else if (it_1->first.second == 1)
            {
                tmp_ui = 2;
            }
            else if (it_1->first.second == 2)
            {
                tmp_ui = 4;
            }
            else
            {
                std::cout << "Undefined Edge NO!\n";
                exit(0);
            }

            oss << JXY[k][JM[k][it_1->first.first](tmp_ui)](0) << ", " << JXY[k][JM[k][it_1->first.first](tmp_ui)](1) << ", " << JXY[k][JM[k][it_1->first.first](tmp_ui)](2) << ";";
            tmp_ui++;
            oss << JXY[k][JM[k][it_1->first.first](tmp_ui)](0) << ", " << JXY[k][JM[k][it_1->first.first](tmp_ui)](1) << ", " << JXY[k][JM[k][it_1->first.first](tmp_ui)](2) << ";";
            tmp_ui++;
            tmp_ui = tmp_ui % (6);
            oss << JXY[k][JM[k][it_1->first.first](tmp_ui)](0) << ", " << JXY[k][JM[k][it_1->first.first](tmp_ui)](1) << ", " << JXY[k][JM[k][it_1->first.first](tmp_ui)](2) << ";";
            it_1++;
        }
        oss << "];\n\n";
        oss << "[mk" << k + 1 << ", nk" << k + 1 << "] = size(JB_2st" << k + 1 << ");\n";
        oss << "for i = 1:mk" << k + 1 << "\n";
        oss << "\txs1 = JB_2st" << k + 1 << "(i, 1)\n";
        oss << "\tys1 = JB_2st" << k + 1 << "(i, 2)\n";
        oss << "\tzs1 = JB_2st" << k + 1 << "(i, 3)\n";
        oss << "\tscatter3(xs1, ys1, zs1, 30, 'k', 'filled')\n";
        oss << "\thold on;\n";
        oss << "end;\n\n";
    }
    //std::cout << "kkk4;\n";
    oss << "\nclear xs1 ys1 zs1;\n\n";

    //< std::vector<std::vector<std::vector<size_t>>> Trace_Node_overall
    for (size_t k = 0; k < Trace_Node_overall.size(); ++k)
    {
        oss << "Tra_Node_" << k + 1 << "=[";
        for (size_t m = 0; m < Trace_Node_overall[k].size(); ++m)
        {
            for (size_t n = 0; n < Trace_Node_overall[k][m].size(); ++n)
            {
                size_t Node_No = Trace_Node_overall[k][m][n];
                oss << JXY[k][Node_No](0) << ", " << JXY[k][Node_No](1) << ", " << JXY[k][Node_No](2) << "; ";
            }
        }
        oss << "];\n";

        oss << "[mo" << k + 1 << ",no" << k + 1 << "] = size(Tra_Node_" << k + 1 << ");\n";
        oss << "for i = 1:mo" << k + 1 << "\n";
        oss << "\txs1 = Tra_Node_" << k + 1 << "(i, 1)\n";
        oss << "\tys1 = Tra_Node_" << k + 1 << "(i, 2)\n";
        oss << "\tzs1 = Tra_Node_" << k + 1 << "(i, 3)\n";
        oss << "\tscatter3(xs1, ys1, zs1, 30, 'k', 'filled')\n";
        oss << "\thold on;\n";
        oss << "end;\n\n";
    }
    //std::cout << "kkk5;\n";
    oss.close();
};

void DFN_mesh::CXX_Triangle_code(string FileKey, const size_t i, const std::vector<Vector3d> tem_verts_trim, const std::vector<std::vector<Vector3d>> Seg_pnt, const double x_ik, const double y_ik, const double element_max)
{
    //Writing data
    string res;
    stringstream ss;
    ss << i;
    ss >> res;

    string B = "_";
    string C = "T";
    FileKey = FileKey + B + res + B + C;
    std::ofstream oss(FileKey, ios::out);

    oss << "Mesh::Unstructured mesh(2);\n";
    int no_trace_seg = 0;
    int no_trace_pnt = 0;
    for (size_t j = 0; j < Seg_pnt.size(); ++j)
    {
        no_trace_pnt += Seg_pnt[j].size();
        no_trace_seg += Seg_pnt[j].size() - 1;
    }
    oss << "mesh.Set(" << tem_verts_trim.size() + no_trace_pnt << ", " << tem_verts_trim.size() + no_trace_seg << ", 1, 0);\n";
    oss << "mesh.SetReg(0, -1, " << element_max << ", " << round(x_ik, 4) << ", " << round(y_ik, 4) << ");\n";

    int domain_pnt_id = 0;
    for (size_t k = 0; k < tem_verts_trim.size(); ++k)
    {

        oss << "mesh.SetPnt(" << k << ", -" << k << ", " << round(tem_verts_trim[k](0), 4) << ", " << round(tem_verts_trim[k](1), 4) << ");"
            << "\n";
        if (k == tem_verts_trim.size() - 1)
            domain_pnt_id = k;
    }
    //set segments, and segment ID must be unique
    int domain_seg_id = 0;
    for (size_t k = 0; k < tem_verts_trim.size(); ++k)
    {
        int k_plus_1 = k + 1 - (int)((k + 1) / tem_verts_trim.size()) * (k + 1);
        oss << "mesh.SetSeg(" << k << ", -" << k << ", " << k << ", " << k_plus_1 << ");\n";
        if (k == tem_verts_trim.size() - 1)
            domain_seg_id = k;
    }

    int intersection_pnt_id = domain_pnt_id + 1;
    int intersection_seg_id = domain_seg_id + 1;
    for (size_t j = 0; j < Seg_pnt.size(); ++j)
    {
        if (Seg_pnt[j].size() >= 2)
        {
            if (Seg_pnt[j].size() == 2 && abs(Seg_pnt[j][0](0) - Seg_pnt[j][1](0)) < 0.0001 && abs(Seg_pnt[j][0](1) - Seg_pnt[j][1](1)) < 0.0001)
            {
                //if the trace is too small, like a point
                oss << "mesh.SetPnt(" << 1 + intersection_pnt_id << ", " << -(1 + intersection_pnt_id) << ", " << round(Seg_pnt[j][0](0), 4) << ", " << round(Seg_pnt[j][0](1), 4) << ");\n";
                intersection_pnt_id++;
            }
            else
            {
                for (size_t k = 0; k < Seg_pnt[j].size(); ++k)
                {
                    oss << "mesh.SetPnt(" << k + intersection_pnt_id << ", " << -(k + intersection_pnt_id) << ", " << round(Seg_pnt[j][k](0), 4) << ", " << round(Seg_pnt[j][k](1), 4) << ");\n";
                    if (k == Seg_pnt[j].size() - 1)
                        intersection_pnt_id += (k + 1);
                };
                for (size_t k = 0; k < Seg_pnt[j].size() - 1; ++k)
                {

                    oss << "mesh.SetSeg(" << k + intersection_seg_id << ", " << -(k + intersection_seg_id) << ", " << k + (intersection_pnt_id - Seg_pnt[j].size()) << ", " << k + (intersection_pnt_id - Seg_pnt[j].size()) + 1 << ");\n";
                    if (k == Seg_pnt[j].size() - 2)
                        intersection_seg_id += (k + 1);
                }
            }
        }
        else
        {
            std::cout << "Error happens in the Class 'Mesh_frac', Function 'DFN_mesh'!\n";
        }
    }

    oss.close();
};

void DFN_mesh::CXX_Fade2D_code(string FileKey, const size_t i, std::map<std::pair<double, double>, int> MapPnt, std::vector<Vector2d> SegIDtoID, const std::vector<Vector3d> tem_verts_trim_AAA)
{
    //std::cout << "---------------Writing data\n";
    //Writing data
    string res;
    stringstream ss;
    ss << i;
    ss >> res;

    string B = "_";
    string C = "F";
    FileKey = FileKey + B + res + B + C;
    std::ofstream oss(FileKey, ios::out);

    /*
        std::map<std::pair<double, double>, int>::iterator it;
        std::map<std::pair<double, double>, int>::iterator itEnd;
        it = MapPnt.begin();
        itEnd = MapPnt.end();
        while (it != itEnd)
        {
            mesh.SetPnt(it->second, -(it->second), (it->first.first), (it->first.second));
            //std::cout << "Pnt ID: " << it->second << ", [" << it->first.first << ", " << it->first.second << "];\n";
            it++;
        }
        //set segments, and segment ID must be unique
        for (int k = 0; k < SegIDtoID.size(); ++k)
        {

            mesh.SetSeg(k, -k, SegIDtoID[k](0), SegIDtoID[k](1));
            //std::cout << "Segment ID: " << k << "; Tag << " << -k << "; Point IDs: " << SegIDtoID[k](0) << ", " << SegIDtoID[k](1) << std::endl;
            //mesh.SetSeg(k, -k, k /point id/, k_plus_1 /point id/);
        }
        */
    oss << "\nstd::vector<Point2> vInputPoints, vInputPoints_A;\n";
    std::map<std::pair<double, double>, int>::iterator it;
    std::map<std::pair<double, double>, int>::iterator itEnd;

    it = MapPnt.begin();
    itEnd = MapPnt.end();
    /*
        double x0, y0;
        int ui = 0, ux = 0;
        for (int k = 0; k < SegIDtoID.size(); ++k)
        {
            if (k != 0 && SegIDtoID[k](1) == 0)
            {
                ui = k;
                ux = SegIDtoID[k](0);
                std::cout << "The last pnt is ID " << ux << std::endl;
            }
        }*/

    oss << "vInputPoints.resize(" << MapPnt.size() << ");\n";

    while (it != itEnd)
    {
        //if (it->second < ux + 1)
        //{
        double x1 = round((it->first.first), 4);
        double y1 = round((it->first.second), 4);
        oss << "vInputPoints[" << it->second << "] = (Point2(" << x1 << ", " << y1 << "));\n";
        //std::cout << "pnt ID = " << it->second << "\n";
        //}
        it++;
    }

    oss << "\nFade_2D dt;\n";
    oss << "dt.insert(vInputPoints);\n";
    oss << "//dt.show(example3_pnt.ps);\n\n";
    oss << "\n";
    /*
        while (it != itEnd)
        {
            double x1 = round((it->first.first), 4);
            double y1 = round((it->first.second), 4);
            oss << "vInputPoints.push_back(Point2(" << x1 << ", " << y1 << "));\n";
            it++;
        }*/
    /*
        std::map<std::pair<double, double>, int>::iterator it_s;
        it_s = MapPnt.begin();
        std::vector<Vector2d> tmp_jxy;
        tmp_jxy.resize(MapPnt.size());
        size_t yu = 0;
        while (it_s != itEnd)
        {
            Vector2d Ax;
            Ax(0) = it_s->first.first;
            Ax(1) = it_s->first.second;
            tmp_jxy[yu](0) = Ax(0);
            tmp_jxy[yu](1) = Ax(1);
            it_s++;
            yu++;
        }*/
    oss << "\nstd::vector<Segment2> vSegments;\n";

    for (size_t k = 0; k < SegIDtoID.size(); ++k)
    {
        oss << "vSegments.push_back(Segment2(vInputPoints[" << SegIDtoID[k](0) << "], vInputPoints[" << SegIDtoID[k](1) << "]));\n";
        //oss << "vSegments.push_back(Segment2(Point2(" << tmp_jxy[SegIDtoID[k](0)](0) << ", " << tmp_jxy[SegIDtoID[k](0)](1) << "),Point2(" << tmp_jxy[SegIDtoID[k](1)](0) << ", " << tmp_jxy[SegIDtoID[k](1)](1) << ")));\n";
    }

    oss.close();
    //std::cout << "---------------finish Writing data\n";
};

void DFN_mesh::Matlab_plot_2D_frac(string FileKey, const size_t i, const std::vector<Vector3d> tem_verts_trim, const std::vector<std::vector<Vector3d>> Seg_pnt)
{
    //Writing data
    string res;
    stringstream ss;
    ss << i;
    ss >> res;

    string B = "_";
    string C = "D";
    FileKey = FileKey + B + res + B + C;
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "verts = [";
    for (size_t i = 0; i < tem_verts_trim.size(); ++i)
    {
        oss << tem_verts_trim[i](0) << ", " << tem_verts_trim[i](1) << "; ";
    }
    oss << "];\n";

    oss << "scatter(verts(:,1), verts(:,2), 30, [rand rand rand]);\n";
    oss << "hold on;\n";

    for (size_t m = 0; m < Seg_pnt.size(); ++m)
    {
        oss << "segment_" << m << " = [";
        for (size_t i = 0; i < Seg_pnt[m].size(); ++i)
        {
            oss << Seg_pnt[m][i](0) << ", " << Seg_pnt[m][i](1) << "; ";
        }
        oss << "];\n";
        oss << "scatter(segment_" << m << "(:,1), segment_" << m << "(:,2), 30, [rand rand rand]);\n";
        oss << "hold on;\n";
    }

    oss.close();
};

inline void DFN_mesh::Address_Repetitive_Node(std::vector<std::vector<Vector3d>> &Coe_Matr_guide,
                                              std::vector<std::vector<Vector3d>> &Coe_Matr_guide_m, const std::vector<std::vector<Vector3d>> JXY, size_t &overall_matrix_dimension)
{
    overall_matrix_dimension = 0;
    NO_Nodes_p = 0;
    Coe_Matr_guide.resize(JXY.size());
    Coe_Matr_guide_m.resize(JXY.size());
    for (size_t i = 0; i < JXY.size(); ++i)
    {
        Coe_Matr_guide[i].resize(JXY[i].size());
        Coe_Matr_guide_m[i].resize(NO_Nodes_p_each_frac[i]);
        for (size_t j = 0; j < JXY[i].size(); ++j)
        {
            int xu, yu;
            bool IY = Find_repetitive_thing(JXY[i][j], JXY, i, j, xu, yu);
            if (IY == true)
            {
                Coe_Matr_guide[i][j](0) = -1;
                Coe_Matr_guide[i][j](1) = xu;
                Coe_Matr_guide[i][j](2) = yu;
            }
            else
            {
                Coe_Matr_guide[i][j](0) = overall_matrix_dimension;
                Coe_Matr_guide[i][j](1) = i;
                Coe_Matr_guide[i][j](2) = j;
                overall_matrix_dimension++;
            }

            if (j < NO_Nodes_p_each_frac[i])
            {
                int xu_S, yu_S;
                bool IY_S = Find_repetitive_thing_p(JXY[i][j], JXY, i, j, xu_S, yu_S, NO_Nodes_p_each_frac);
                if (IY_S == true)
                {
                    Coe_Matr_guide_m[i][j](0) = -1;
                    Coe_Matr_guide_m[i][j](1) = xu_S;
                    Coe_Matr_guide_m[i][j](2) = yu_S;
                }
                else
                {
                    Coe_Matr_guide_m[i][j](0) = NO_Nodes_p;
                    Coe_Matr_guide_m[i][j](1) = i;
                    Coe_Matr_guide_m[i][j](2) = j;
                    NO_Nodes_p++;
                }
            }
            //std::cout << JXY[i][j] << ", " << IY <<"; ";
        }
        //std::cout <<std::endl;
    }
};

inline void DFN_mesh::Remove_unnecessary_fractures(size_t Cluster_ID, DFN::Domain dom)
{
    //remove fractures that have no contribution to flow

    //--------test
    /*
    for (size_t i = 0; i < Listofclusters_mesh_only[Cluster_ID].size(); ++i)
    {
        int frac_ID = Listofclusters_mesh_only[Cluster_ID][i];
        if (dom.Fractures[frac_ID].If_intersect_surfaces[0] == 1)
        {
            std::cout << "Fracture ID = " << frac_ID << " intersects domain frace 0 (top)\n";
        };
    }

    for (size_t i = 0; i < Listofclusters_mesh_only[Cluster_ID].size(); ++i)
    {
        int frac_ID = Listofclusters_mesh_only[Cluster_ID][i];
        if (dom.Fractures[frac_ID].If_intersect_surfaces[1] == 1)
        {
            std::cout << "Fracture ID = " << frac_ID << " intersects domain frace 1 (bottom)\n";
        };
    }*/
    //----------
    /*
    std::cout << "Percolation cluster: ";
    for (size_t i = 0; i < Listofclusters_mesh_only[Cluster_ID].size(); ++i)
    {
        std::cout << Listofclusters_mesh_only[Cluster_ID][i] << ", ";
    }
    std::cout << "\n";*/

    size_t i = 0;
    std::vector<size_t> BC_tag = {1, 1, 0, 0, 0, 0};
    while (i < Listofclusters_mesh_only[Cluster_ID].size())
    {
        int fracID = Listofclusters_mesh_only[Cluster_ID][i];
        //std::cout << "Frac (ID = " << fracID << ") is being checking\n";
        if (fracID != -1)
        {
            size_t Frac_in_c = 0;

            for (size_t j = 0; j < BC_tag.size(); ++j)
            {
                if (BC_tag[j] == 1)
                {
                    if (dom.Fractures[fracID].If_intersect_surfaces[j] == 1)
                    {
                        // std::cout << "\t//It intersects domain Face " << j << ";\n";
                        Frac_in_c++;
                    }
                }
            }

            for (size_t j = 0; j < Listofclusters_mesh_only[Cluster_ID].size(); ++j)
            {
                if (j != i && Listofclusters_mesh_only[Cluster_ID][j] >= 0)
                {
                    // std::set<size_t> Intersect_other_frac_after_trim;
                    /*
                    std::set<size_t>::iterator it;
                    it = dom.Fractures[Listofclusters_mesh_only[Cluster_ID][j]].Intersect_other_frac_after_trim.begin();
                    std::cout << "\tFracture ID " << Listofclusters_mesh_only[Cluster_ID][j] << " intersect with ";
                    while (it != dom.Fractures[Listofclusters_mesh_only[Cluster_ID][j]].Intersect_other_frac_after_trim.end())
                    {
                        std::cout << *it << ", ";
                        it++;
                    }
                    std::cout << "\n";
                    */
                    std::set<size_t>::iterator pp = dom.Fractures[Listofclusters_mesh_only[Cluster_ID][j]].Intersect_other_frac_after_trim.find(fracID);
                    if (pp != dom.Fractures[Listofclusters_mesh_only[Cluster_ID][j]].Intersect_other_frac_after_trim.end())
                    {
                        //std::cout << "\t**This fracture (ID = " << fracID << ") intersects Fracture ID " << Listofclusters_mesh_only[Cluster_ID][j] << ";\n";
                        Frac_in_c++;
                    }
                }
            }

            if (Frac_in_c < 2)
            {
                //std::cout << "\tremove fracture ID: " << Listofclusters_mesh_only[Cluster_ID][i] << std::endl;
                //std::cout << std::endl;
                Listofclusters_mesh_only[Cluster_ID][i] = -1;
                i = 0;
                //std::cout << "Now checking again!\n";
            }
            else
            {
                ++i;
            }
        }
        else
        {
            ++i;
            //std::cout << "\tit is a removed fracture!\n";
        }
    }
    //----------------
};

} // namespace DFN