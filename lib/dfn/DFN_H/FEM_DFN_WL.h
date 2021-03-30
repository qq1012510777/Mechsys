#pragma once
#include "../Math_WL_H/Math_WL.h"
#include "Dense"
#include "Mesh_DFN_WL.h"
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{
class FEM_DFN
{
public:
    //std::vector<std::vector<std::pair<Vector3d, Vector3d>>> Flux_rate;
    //std::vector<std::vector<std::pair<Vector3d, Vector3d>>> Flux_rate_2D;
    ///<store the coordinates (x,y,z) that represent an element (usually it is the center of the element), and the flux rate vector at this point

    std::vector<std::vector<std::vector<std::pair<std::pair<Vector3d, Vector3d>, std::pair<Vector3d, Vector3d>>>>> Flux_rate_Bound_in;
    std::vector<std::vector<std::vector<std::pair<std::pair<Vector3d, Vector3d>, std::pair<Vector3d, Vector3d>>>>> Flux_rate_Bound_out;
    ///< 2d and 3d, (pnt, q)

    //std::vector<std::vector<RowVector6d>> pd_N_pd_x;
    // also the gridient

    //std::vector<std::vector<RowVector6d>> pd_N_pd_y;

    std::vector<std::vector<Vector3d>> Velo_3D;

public:
    FEM_DFN(DFN::DFN_mesh DFN_mesh,
            DFN::Domain dom);

    void Assemble_matrix(DFN::DFN_mesh DFN_mesh,
                         const double Kper,
                         double *K_overall,
                         double *F_overall);

    void Matlab_FEM_head_result(string FileKey,
                                DFN::DFN_mesh DFN_mesh,
                                double *X_overall,
                                DFN::Domain dom);

    void Verify_in_and_out_flux(DFN::DFN_mesh DFN_mesh,
                                double *X_overall);

    void Calculate_flux_of_A_Boun_ele(DFN::DFN_mesh DFN_mesh,
                                      double *X_overall,
                                      double K_coe,
                                      size_t s_ifrac,
                                      size_t s_jele,
                                      size_t EdgeNO,
                                      std::vector<std::pair<std::pair<Vector3d, Vector3d>, std::pair<Vector3d, Vector3d>>> &SF,
                                      size_t NumOfPnt_q);

    void Weighted_flux_shared_pnt(DFN::DFN_mesh DFN_mesh,
                                  double *X_overall,
                                  double K_coe,
                                  double A_ele,
                                  double &qx,
                                  double &qy,
                                  size_t ty_ifrac,
                                  size_t ty_jele,
                                  size_t ty_edge_no,     // 0 - 2
                                  size_t ty_local_pnt_no // 0 - 2
    );

    void Recursion_find_ele_shared_a_pnt(DFN::DFN_mesh DFN_mesh,
                                         size_t ifrac,
                                         size_t jele,
                                         size_t pnt_ID_p,
                                         size_t local_pnt_ID_p,
                                         std::map<size_t, Vector2d> &AAD);

    void Calculated_velocity_of_a_pnt_local(double xi,
                                            double eta,
                                            double K_coe,
                                            VectorXd h_e,
                                            VectorXd x_e,
                                            VectorXd y_e,
                                            double &qx,
                                            double &qy);

    void Calculated_pnt_from_local_to_golbal(double xi,
                                             double eta,
                                             VectorXd x_e,
                                             VectorXd y_e,
                                             Vector3d &Center);

    void p_PHI_over_p_x_and_y(double xi,
                              double eta,
                              MatrixXd JXYe_x,
                              MatrixXd JXYe_y,
                              VectorXd &pd_N_over_pd_x,
                              VectorXd &pd_N_over_pd_y,
                              MatrixXd &Jacobi);

    void PHI_shape_function(double xi,
                            double eta,
                            double &N_1,
                            double &N_2,
                            double &N_3,
                            double &N_4,
                            double &N_5,
                            double &N_6);

    void p_PHI_over_p_xi_and_eta(double xi,
                                 double eta,
                                 double &N_natural_a_1_s,
                                 double &N_natural_a_2_s,
                                 double &N_natural_a_3_s,
                                 double &N_natural_a_4_s,
                                 double &N_natural_a_5_s,
                                 double &N_natural_a_6_s,
                                 double &N_natural_b_1_s,
                                 double &N_natural_b_2_s,
                                 double &N_natural_b_3_s,
                                 double &N_natural_b_4_s,
                                 double &N_natural_b_5_s,
                                 double &N_natural_b_6_s);

    void p_PSI_over_p_x_and_y(double xi,
                              double eta,
                              MatrixXd JXYe_x,
                              MatrixXd JXYe_y,
                              VectorXd &pd_N_over_pd_x,
                              VectorXd &pd_N_over_pd_y,
                              MatrixXd &Jacobi);

    void PSI_shape_function(double xi,
                            double eta,
                            double &N_1,
                            double &N_2,
                            double &N_3);

    void p_PSI_over_p_xi_and_eta(double xi,
                                 double eta,
                                 double &N_natural_a_1_s,
                                 double &N_natural_a_2_s,
                                 double &N_natural_a_3_s,
                                 double &N_natural_b_1_s,
                                 double &N_natural_b_2_s,
                                 double &N_natural_b_3_s);

    void Three_D_velocity(DFN::DFN_mesh DFN_mesh, double *X_overall);
};

//*******************************************
inline FEM_DFN::FEM_DFN(DFN::DFN_mesh DFN_mesh, DFN::Domain dom)
{
    double k_coe = 1.0; // permeability of each fracture, relating to aperture

    size_t NUM_NODES_velocity = DFN_mesh.NO_all_pnts;
    size_t NUM_NODES_p = DFN_mesh.NO_Nodes_p;
    size_t Matrix_D = (NUM_NODES_velocity * 2 + NUM_NODES_p);

    double *K_overall = new double[Matrix_D * Matrix_D];
    if (K_overall == NULL)
    {
        cout << "Error! Cannot alloc to matrix 'K_overall'!\n";
        exit(0);
    }
    else
    {
        for (size_t i = 0; i < Matrix_D * Matrix_D; ++i)
        {
            K_overall[i] = 0;
        }
    }

    double *F_overall = new double[Matrix_D];
    if (F_overall == NULL)
    {
        cout << "Error! Cannot alloc to matrix 'F_overall'!\n";
        exit(0);
    }
    else
    {
        for (size_t i = 0; i < Matrix_D; ++i)
        {
            F_overall[i] = 0;
        }
    }

    Assemble_matrix(DFN_mesh, k_coe, K_overall, F_overall);
    /*
    cout << "\nK_overall;\n";
    for (size_t i = 0; i < Matrix_D; ++i)
    {
        for (size_t j = 0; j < Matrix_D; ++j)
        {
            size_t Idx = i * Matrix_D + j;

            cout << K_overall[Idx] << ",";
        }
        cout << "|" << endl;
    }

    cout << "\nF_overall;\n";
    for (size_t i = 0; i < Matrix_D; ++i)
    {
        cout << F_overall[i] << endl;
    }*/

    // K x = F
    int n = Matrix_D; // dimensions of coefficient matrix K (2D)
    int m = 1;        // dimensions of F matrix
    //double a = 1, b = -1.00;

    //dgemv_("N", &n, &n, &a, K_overall, &n, F_overall, &m, &b, X_overall, &m);
    int *ipiv = new int[n];
    if (ipiv == NULL)
    {
        cout << "Error! Cannot alloc to matrix 'F_overall'!\n";
        exit(0);
    }
    else
    {
        for (size_t i = 0; i < (size_t)n; ++i)
        {
            ipiv[i] = 0;
        }
    }
    int info;
    std::cout << "start solving matrix;\n";
    dgesv_(&n, &m, K_overall, &n, ipiv, F_overall, &n, &info);
    std::cout << "finish solving matrix;\n";
    // now F_overall is the solution X
    if (info > 0)
    {
        cout << "The solution could not be computed!!\n";
        exit(0);
    }

    /*cout << "\nn " << n << endl;
    cout << "a " << a << endl;
    cout << "b " << b << endl;*/

    /*
    cout << "\nX_overall with BLAS;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        cout << X_overall[i] << endl;
    }
    */

    /*
    cout << "\n22 K_overall;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        for (size_t j = 0; j < DFN_mesh.Matrix_dimesions; ++j)
        {
            size_t Idx = i * DFN_mesh.Matrix_dimesions + j;

            cout << K_overall[Idx] << "\t";
            size_t Idy = j * DFN_mesh.Matrix_dimesions + i;
            if (abs(K_overall[Idx] - K_overall[Idy]) > 0.0001)
            {
                cout << "Error! \n";
                exit(0);
            }
        }
        cout << "|" << endl;
    }

    cout << "\n22 F_overall;\n";
    for (size_t i = 0; i < DFN_mesh.Matrix_dimesions; ++i)
    {
        cout << F_overall[i] << endl;
    }
    */
    //dgemv_("N", DFN_mesh.Matrix_dimesions, DFN_mesh.Matrix_dimesions, 1, K_overall, /*LDA*/ DFN_mesh.Matrix_dimesions, X_overall, DFN_mesh.Matrix_dimesions /*INCX*/, 1, F_overall, DFN_mesh.Matrix_dimesions);

    Verify_in_and_out_flux(DFN_mesh, F_overall);
    Three_D_velocity(DFN_mesh, F_overall);
    Matlab_FEM_head_result("tdfn_color_head.m", DFN_mesh, F_overall, dom);
    delete[] ipiv;
    ipiv = NULL;
    delete[] F_overall;
    F_overall = NULL;
    delete[] K_overall;
    K_overall = NULL;
};
//*******************************************

inline void FEM_DFN::Assemble_matrix(DFN::DFN_mesh DFN_mesh, const double Kper, double *K_overall, double *F_overall)
{
    size_t NUM_NODES_velocity = DFN_mesh.NO_all_pnts;
    size_t NUM_NODES_p = DFN_mesh.NO_Nodes_p;
    size_t Matrix_D = (NUM_NODES_velocity * 2 + NUM_NODES_p);
    /*
    Eigen::MatrixXd K, F;
    K = Eigen::MatrixXd::Zero(DFN_mesh.Matrix_dimesions, DFN_mesh.Matrix_dimesions);
    F = Eigen::MatrixXd::Zero(DFN_mesh.Matrix_dimesions, 1);*/
    //pd_N_pd_x.resize(DFN_mesh.JM.size());
    //pd_N_pd_y.resize(DFN_mesh.JM.size());
    for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
    {
        //pd_N_pd_x[i].resize(DFN_mesh.JM[i].size()); // the fracture NO. i
        //pd_N_pd_y[i].resize(DFN_mesh.JM[i].size());
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            Eigen::RowVectorXd w, xi, eta;
            w = Eigen::RowVectorXd::Zero(6);
            w << 0.1713244923791700, 0.3607615730481380, 0.4679139345726910, 0.1713244923791700, 0.3607615730481380, 0.4679139345726910;
            xi = Eigen::RowVectorXd::Zero(6);
            xi << 0.9324695142031520, 0.6612093864662640, 0.2386191860831960, -0.9324695142031520, -0.6612093864662640, -0.2386191860831960;
            eta = xi;

            MatrixXd D1e;
            D1e = MatrixXd::Zero(6, 6);

            for (size_t ik = 0; ik < (size_t)D1e.rows(); ++ik)
            {
                for (size_t jk = 0; jk < (size_t)D1e.cols(); ++jk)
                {
                    Eigen::MatrixXd JXYe_x, JXYe_y;
                    JXYe_x = Eigen::MatrixXd::Zero(6, 1);
                    JXYe_y = Eigen::MatrixXd::Zero(6, 1);

                    for (size_t kk = 0; kk < (size_t)JXYe_x.rows(); ++kk)
                    {
                        JXYe_x(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](kk)](0);
                        JXYe_y(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](kk)](1);
                    }

                    VectorXd pd_N_over_pd_x, pd_N_over_pd_y;
                    pd_N_over_pd_x = VectorXd::Zero(6);
                    pd_N_over_pd_y = VectorXd::Zero(6);

                    MatrixXd Jacobi;
                    Jacobi = MatrixXd::Zero(2, 2);

                    Eigen::VectorXd PHI;
                    PHI = Eigen::VectorXd::Zero(6);

                    double Phi_1;
                    double Phi_2;
                    double Phi_3;
                    double Phi_4;
                    double Phi_5;
                    double Phi_6;
                    PHI_shape_function(xi(ik),
                                       eta(jk),
                                       Phi_1,
                                       Phi_2,
                                       Phi_3,
                                       Phi_4,
                                       Phi_5,
                                       Phi_6);
                    PHI << Phi_1, Phi_4, Phi_2, Phi_5, Phi_3, Phi_6;

                    p_PHI_over_p_x_and_y(xi(ik), eta(jk), JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

                    D1e = D1e + 1. / Kper * w(ik) * w(jk) * Jacobi.determinant() * PHI * PHI.transpose();
                };
            };

            MatrixXd C1e, C2e;
            C1e = MatrixXd::Zero(6, 3);
            C2e = MatrixXd::Zero(6, 3);
            for (size_t ik = 0; ik < 6; ++ik)
            {
                for (size_t jk = 0; jk < 6; ++jk)
                {
                    Eigen::VectorXd PHI;
                    PHI = Eigen::VectorXd::Zero(6);

                    double Phi_1;
                    double Phi_2;
                    double Phi_3;
                    double Phi_4;
                    double Phi_5;
                    double Phi_6;
                    PHI_shape_function(xi(ik),
                                       eta(jk),
                                       Phi_1,
                                       Phi_2,
                                       Phi_3,
                                       Phi_4,
                                       Phi_5,
                                       Phi_6);
                    PHI << Phi_1, Phi_4, Phi_2, Phi_5, Phi_3, Phi_6;
                    //cout << "0.021\n";
                    Eigen::MatrixXd JXYe_x, JXYe_y;
                    JXYe_x = Eigen::MatrixXd::Zero(3, 1);
                    JXYe_y = Eigen::MatrixXd::Zero(3, 1);
                    //cout << "0.0214\n";
                    for (size_t kk = 0; kk < 3; kk++)
                    {
                        int ff = -1;
                        if (kk == 0)
                            ff = 0;
                        else if (kk == 1)
                            ff = 2;
                        else if (kk == 2)
                            ff = 4;

                        JXYe_x(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](ff)](0);
                        JXYe_y(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](ff)](1);
                    }
                    //cout << "0.0215\n";
                    VectorXd pd_N_over_pd_x, pd_N_over_pd_y;
                    pd_N_over_pd_x = VectorXd::Zero(3);
                    pd_N_over_pd_y = VectorXd::Zero(3);

                    MatrixXd Jacobi;
                    Jacobi = MatrixXd::Zero(2, 2);
                    //cout << "0.022\n";

                    p_PSI_over_p_x_and_y(xi(ik), eta(jk), JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

                    C1e = C1e + w(ik) * w(jk) * PHI * pd_N_over_pd_x.transpose() * Jacobi.determinant();
                    C2e = C2e + w(ik) * w(jk) * PHI * pd_N_over_pd_y.transpose() * Jacobi.determinant();
                }
            }

            MatrixXd B1e, B2e;
            B1e = MatrixXd::Zero(3, 6);
            B2e = MatrixXd::Zero(3, 6);
            for (size_t ik = 0; ik < 6; ++ik)
            {
                for (size_t jk = 0; jk < 6; ++jk)
                {
                    Eigen::MatrixXd JXYe_x, JXYe_y;
                    JXYe_x = Eigen::MatrixXd::Zero(3, 1);
                    JXYe_y = Eigen::MatrixXd::Zero(3, 1);

                    for (size_t kk = 0; kk < 3; kk++)
                    {
                        int ff = -1;
                        if (kk == 0)
                            ff = 0;
                        else if (kk == 1)
                            ff = 2;
                        else if (kk == 2)
                            ff = 4;

                        JXYe_x(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](ff)](0);
                        JXYe_y(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](ff)](1);
                    }

                    VectorXd pd_N_over_pd_x, pd_N_over_pd_y;
                    pd_N_over_pd_x = VectorXd::Zero(3);
                    pd_N_over_pd_y = VectorXd::Zero(3);

                    MatrixXd Jacobi;
                    Jacobi = MatrixXd::Zero(2, 2);

                    p_PSI_over_p_x_and_y(xi(ik), eta(jk), JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

                    Eigen::VectorXd PHI;
                    PHI = Eigen::VectorXd::Zero(6);

                    double Phi_1;
                    double Phi_2;
                    double Phi_3;
                    double Phi_4;
                    double Phi_5;
                    double Phi_6;
                    PHI_shape_function(xi(ik),
                                       eta(jk),
                                       Phi_1,
                                       Phi_2,
                                       Phi_3,
                                       Phi_4,
                                       Phi_5,
                                       Phi_6);
                    PHI << Phi_1, Phi_4, Phi_2, Phi_5, Phi_3, Phi_6;

                    B1e = B1e + w(ik) * w(jk) * pd_N_over_pd_x * PHI.transpose() * Jacobi.determinant();
                    B2e = B2e + w(ik) * w(jk) * pd_N_over_pd_y * PHI.transpose() * Jacobi.determinant();
                }
            }

            //-------now assemble the K matrix
            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    size_t Node_m = DFN_mesh.JM[i][j](jq);
                    size_t Node_n = DFN_mesh.JM[i][j](kq);
                    if (DFN_mesh.Coe_Matr_guide[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Coe_Matr_guide[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](kq)](0);
                    }

                    size_t Global_Idx = Node_m * Matrix_D + Node_n;

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(D1e(jq, kq), 4);
                }
            }

            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    size_t Node_m = DFN_mesh.JM[i][j](jq);
                    size_t Node_n = DFN_mesh.JM[i][j](kq);
                    if (DFN_mesh.Coe_Matr_guide[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Coe_Matr_guide[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](kq)](0);
                    }

                    size_t Global_Idx = (Node_m + NUM_NODES_velocity) * Matrix_D + (Node_n + NUM_NODES_velocity);

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(D1e(jq, kq), 4);
                }
            }

            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 3; ++kq)
                {
                    int yt = -1;
                    if (kq == 0)
                        yt = 0;
                    else if (kq == 1)
                        yt = 2;
                    else if (kq == 2)
                        yt = 4;

                    size_t Node_m = DFN_mesh.JM[i][j](jq);
                    size_t Node_n = DFN_mesh.JM[i][j](yt);
                    if (DFN_mesh.Coe_Matr_guide[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Coe_Matr_guide_m[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Coe_Matr_guide_m[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide_m[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Coe_Matr_guide_m[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Coe_Matr_guide_m[i][DFN_mesh.JM[i][j](yt)](0);
                    }

                    size_t Global_Idx = (Node_m)*Matrix_D + (Node_n + 2 * NUM_NODES_velocity);

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(C1e(jq, kq), 4);
                }
            }

            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 3; ++kq)
                {
                    int yt = -1;
                    if (kq == 0)
                        yt = 0;
                    else if (kq == 1)
                        yt = 2;
                    else if (kq == 2)
                        yt = 4;

                    size_t Node_m = DFN_mesh.JM[i][j](jq);
                    size_t Node_n = DFN_mesh.JM[i][j](yt);
                    if (DFN_mesh.Coe_Matr_guide[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Coe_Matr_guide_m[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Coe_Matr_guide_m[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide_m[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Coe_Matr_guide_m[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Coe_Matr_guide_m[i][DFN_mesh.JM[i][j](yt)](0);
                    }

                    size_t Global_Idx = (Node_m + NUM_NODES_velocity) * Matrix_D + (Node_n + 2 * NUM_NODES_velocity);

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(C2e(jq, kq), 4);
                }
            }

            for (size_t jq = 0; jq < 3; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    int yt = -1;
                    if (jq == 0)
                        yt = 0;
                    else if (jq == 1)
                        yt = 2;
                    else if (jq == 2)
                        yt = 4;

                    size_t Node_m = DFN_mesh.JM[i][j](yt);
                    size_t Node_n = DFN_mesh.JM[i][j](kq);
                    if (DFN_mesh.Coe_Matr_guide_m[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Coe_Matr_guide_m[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide_m[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Coe_Matr_guide_m[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Coe_Matr_guide_m[i][DFN_mesh.JM[i][j](yt)](0);
                    }

                    if (DFN_mesh.Coe_Matr_guide[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](kq)](0);
                    }

                    size_t Global_Idx = (Node_m + 2 * NUM_NODES_velocity) * Matrix_D + (Node_n);

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(B1e(jq, kq), 4);
                }
            }

            for (size_t jq = 0; jq < 3; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    int yt = -1;
                    if (jq == 0)
                        yt = 0;
                    else if (jq == 1)
                        yt = 2;
                    else if (jq == 2)
                        yt = 4;

                    size_t Node_m = DFN_mesh.JM[i][j](yt);
                    size_t Node_n = DFN_mesh.JM[i][j](kq);
                    if (DFN_mesh.Coe_Matr_guide_m[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Coe_Matr_guide_m[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide_m[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Coe_Matr_guide_m[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Coe_Matr_guide_m[i][DFN_mesh.JM[i][j](yt)](0);
                    }

                    if (DFN_mesh.Coe_Matr_guide[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Coe_Matr_guide[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Coe_Matr_guide[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Coe_Matr_guide[i][DFN_mesh.JM[i][j](kq)](0);
                    }

                    size_t Global_Idx = (Node_m + 2 * NUM_NODES_velocity) * Matrix_D + (Node_n + DFN_mesh.NO_all_pnts);

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(B2e(jq, kq), 4);
                }
            }
        }

        std::map<std::pair<size_t, size_t>, Vector5d>::iterator it;
        it = DFN_mesh.JB_2[i].begin();
        while (it != DFN_mesh.JB_2[i].end())
        {
            //std::vector<std::map<std::pair<size_t, size_t>, Vector5d>> DFN_mesh.JB_2;
            int II = it->first.first;                   //DFN_mesh.JB_2[i](ti, 0);
            int bondary_side_number = it->first.second; //DFN_mesh.JB_2[i](ti, 1);

            //Eigen::VectorXd Fe;
            //Fe = Eigen::VectorXd::Zero(6);

            double G1, G2;
            Eigen::RowVectorXd Jx, Jy;
            Jx = Eigen::RowVectorXd::Zero(6);
            Jy = Jx;
            for (size_t j = 0; j < 6; ++j)
            {
                Jx(j) = DFN_mesh.JXY[i][DFN_mesh.JM[i][II](j)](0);
                Jy(j) = DFN_mesh.JXY[i][DFN_mesh.JM[i][II](j)](1);
            }

            double q0 = it->second(2);
            double q1 = it->second(3);
            double q2 = it->second(4);

            if (bondary_side_number == 0)
            {
                double L = pow(pow(Jx(0) - Jx(2), 2) + pow(Jy(0) - Jy(2), 2), 0.5);
                L = L * 0.5;

                G1 = (L * (q0 + 2 * q1)) / 3;
                G2 = (L * (2 * q1 + q2)) / 3;
            }
            else if (bondary_side_number == 1)
            {
                double L = pow(pow(Jx(2) - Jx(4), 2) + pow(Jy(2) - Jy(4), 2), 0.5);
                L = L * 0.5;

                G1 = (L * (q0 + 2 * q1)) / 3;
                G2 = (L * (2 * q1 + q2)) / 3;
            }
            else if (bondary_side_number == 2)
            {
                double L = pow(pow(Jx(4) - Jx(0), 2) + pow(Jy(4) - Jy(0), 2), 0.5);
                L = L * 0.5;

                G1 = (L * (q0 + 2 * q1)) / 3;
                G2 = (L * (2 * q1 + q2)) / 3;
            }

            Vector2d Ge;
            Ge << round(G1, 4), round(G2, 4);
            //cout << "Fe: " << ti << endl;

            for (size_t s = 0, si = 0; s < 3; s = s + 2, si++)
            {
                size_t Node_s = DFN_mesh.JM[i][II]((bondary_side_number * 2 + s) % 6);
                if (DFN_mesh.Coe_Matr_guide_m[i][Node_s](0) == -1)
                {
                    size_t i_frac = DFN_mesh.Coe_Matr_guide_m[i][Node_s](1); // the i th frac
                    size_t j_node = DFN_mesh.Coe_Matr_guide_m[i][Node_s](2); // the j th node
                    Node_s = DFN_mesh.Coe_Matr_guide_m[i_frac][j_node](0);
                }
                else
                {
                    Node_s = DFN_mesh.Coe_Matr_guide_m[i][DFN_mesh.JM[i][II](s)](0);
                }
                //cout << Fe(s, 0) << endl;

                F_overall[Node_s + NUM_NODES_velocity * 2] = F_overall[Node_s + NUM_NODES_velocity * 2] + Ge(si);
                //F(JM(II - 1, s) - 1, 0) = F(JM(II - 1, s) - 1, 0) + Fe(s, 0);
                //F(Node_s, 0) = F(Node_s, 0) + Fe(s, 0);
            }

            it++;
        }

        std::map<size_t, double>::iterator it_s;
        it_s = DFN_mesh.JB_1[i].begin();
        while (it_s != DFN_mesh.JB_1[i].end())
        {
            size_t NODE_s = it_s->first;
            double BC_1 = it_s->second; // node pressure value

            size_t IDX; // node ID
            int Guide_1 = DFN_mesh.Coe_Matr_guide_m[i][NODE_s](0);
            if (Guide_1 == -1)
            {
                size_t frac_i = DFN_mesh.Coe_Matr_guide_m[i][NODE_s](1);
                size_t node_j = DFN_mesh.Coe_Matr_guide_m[i][NODE_s](2);
                IDX = DFN_mesh.Coe_Matr_guide_m[frac_i][node_j](0);
            }
            else
            {
                IDX = DFN_mesh.Coe_Matr_guide_m[i][NODE_s](0);
            }

            for (size_t go = 0; go < Matrix_D; ++go)
            {
                size_t golbalID = go * Matrix_D + IDX + 2 * NUM_NODES_velocity;
                F_overall[go] = F_overall[go] - K_overall[golbalID] * BC_1;

                K_overall[golbalID] = 0;

                size_t golbalID_row = (IDX + 2 * NUM_NODES_velocity) * Matrix_D + go;
                K_overall[golbalID_row] = 0;
            }

            K_overall[(IDX + 2 * NUM_NODES_velocity) * Matrix_D + (IDX + 2 * NUM_NODES_velocity)] = 1;
            F_overall[(IDX + 2 * NUM_NODES_velocity)] = BC_1;
            it_s++;
        }
    };
};

inline void FEM_DFN::Matlab_FEM_head_result(string FileKey, DFN::DFN_mesh DFN_mesh, double *X_overall, DFN::Domain dom)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n\n";

    size_t NO_frac = DFN_mesh.JXY.size();
    for (size_t i = 0; i < NO_frac; ++i)
    {
        oss << "JXY_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JXY_3D[i].size(); ++j)
        {
            oss << round(DFN_mesh.JXY_3D[i][j](0), 3) << ", ";
            oss << round(DFN_mesh.JXY_3D[i][j](1), 3) << ", ";
            oss << round(DFN_mesh.JXY_3D[i][j](2), 3) << "; ";
        }
        oss << "];\n";

        oss << "JXY_2D_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JXY_3D[i].size(); ++j)
        {
            oss << round(DFN_mesh.JXY[i][j](0), 3) << ", ";
            oss << round(DFN_mesh.JXY[i][j](1), 3) << ", ";
            oss << round(DFN_mesh.JXY[i][j](2), 3) << "; ";
        }
        oss << "];\n";

        oss << "JM_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            oss << DFN_mesh.JM[i][j](0) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](1) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](2) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](3) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](4) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](5) + 1 << "; ";
        }
        oss << "];\n";

        oss << "JXYs_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.NO_Nodes_p_each_frac[i]; ++j)
        {
            oss << round(DFN_mesh.JXY_3D[i][j](0), 3) << ", ";
            oss << round(DFN_mesh.JXY_3D[i][j](1), 3) << ", ";
            oss << round(DFN_mesh.JXY_3D[i][j](2), 3) << "; ";
        }
        oss << "];\n";

        oss << "JXY_2Ds_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.NO_Nodes_p_each_frac[i]; ++j)
        {
            oss << round(DFN_mesh.JXY[i][j](0), 3) << ", ";
            oss << round(DFN_mesh.JXY[i][j](1), 3) << ", ";
            oss << round(DFN_mesh.JXY[i][j](2), 3) << "; ";
        }
        oss << "];\n";

        oss << "JMs_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            oss << DFN_mesh.JM[i][j](0) + 1 << ", ";
            //oss << DFN_mesh.JM[i][j](1) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](2) + 1 << ", ";
            //oss << DFN_mesh.JM[i][j](3) + 1 << ", ";
            oss << DFN_mesh.JM[i][j](4) + 1 << "; ";
            //oss << DFN_mesh.JM[i][j](5) + 1 << "; ";
        }
        oss << "];\n";

        oss << "Data_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.Coe_Matr_guide_m[i].size(); ++j)
        {
            size_t Node_s = 0;
            if (DFN_mesh.Coe_Matr_guide_m[i][j](0) == -1)
            {
                size_t i_frac = DFN_mesh.Coe_Matr_guide_m[i][j](1);
                size_t j_node = DFN_mesh.Coe_Matr_guide_m[i][j](2);
                Node_s = DFN_mesh.Coe_Matr_guide_m[i_frac][j_node](0);
            }
            else
            {
                Node_s = DFN_mesh.Coe_Matr_guide_m[i][j](0);
            }
            oss << round(X_overall[Node_s + 2 * DFN_mesh.NO_all_pnts], 4) << "; ";
        }
        oss << "];\n";
        oss << "figure(1)\n";
        oss << "P_" << i + 1 << " = patch('Vertices', JXYs_" << i + 1 << ", 'Faces', JMs_" << i + 1 << ", 'FaceVertexCData', Data_" << i + 1 << ", 'FaceColor', 'interp', 'EdgeAlpha', 1);\n";
        oss << "hold on;\n";
    }
    oss << "view(3);\n";
    oss << "colorbar;\n";

    oss << "hold on;\n";
    //Plotting the model domain
    oss << "figure(1)\n";
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
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN head distribution');\nhold on;\n";

    oss << endl;
    oss << "%%*****flux rate vector***\n";
    // std::vector<std::vector<std::pair<Vector3d, Vector3d>>> Flux_rate;
    oss << "figure(1)\n";
    for (size_t i = 0; i < NO_frac; ++i)
    {
        oss << "q_vector_" << i + 1 << "=[";
        for (size_t j = 0; j < DFN_mesh.JXY_3D[i].size(); ++j)
        {
            for (size_t yj = 0; yj < 3; ++yj)
                oss << DFN_mesh.JXY_3D[i][j](yj) << ", ";

            for (size_t yj = 0; yj < 3; ++yj)
            {
                //-----
                oss << Velo_3D[i][j](yj);
                if (yj == 2)
                    oss << "; ";
                else
                    oss << ", ";
            }
        }
        oss << "];\n";
        //oss << "quiver3(q_vector_" << i + 1 << "(:,1), q_vector_" << i + 1 << "(:,2), q_vector_" << i + 1 << "(:,3), q_vector_" << i + 1 << "(:,4), q_vector_" << i + 1 << "(:,5), q_vector_" << i + 1 << "(:,6), 'AutoScale','off');\n\n";
        //oss << "hold on;\n";
    };

    oss << "q_vector_3d_all = [";
    for (size_t i = 0; i < NO_frac; ++i)
    {
        oss << "q_vector_" << i + 1 << "; ";
    }
    oss << "];\n";
    oss << "quiver3(q_vector_3d_all(:,1), q_vector_3d_all(:,2), q_vector_3d_all(:,3), q_vector_3d_all(:,4), q_vector_3d_all(:,5), q_vector_3d_all(:,6));\n\n";
    oss << "hold on;\n";

    oss << "\n\n\n%show frac pressure contour and flow rate vectors in 2D\n";

    for (size_t i = 0; i < NO_frac; ++i)
    {
        size_t yk = 0;
        std::map<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>>::iterator
            AS;
        for (size_t f = 0; f < DFN_mesh.Inlet.size(); ++f)
        {
            AS = DFN_mesh.Inlet[f].begin();
            if (AS->first.first == i)
            {
                yk = 1;
            }
        }
        if (yk == 1)
        {
            oss << "\n\n";
            oss << "figure(" << i + 2 << ")\n";

            oss << "x_tr" << i + 1 << " = JXY_2Ds_" << i + 1 << "(:, 1);\n";
            oss << "y_tr" << i + 1 << " = JXY_2Ds_" << i + 1 << "(:, 2);\n";
            oss << "z_tr" << i + 1 << " = Data_" << i + 1 << "(:);\n";

            oss << "nx" << i + 1 << " = 500;\n";
            oss << "ny" << i + 1 << " = 500;\n";
            oss << "[X" << i + 1 << ",Y" << i + 1 << "]";
            oss << " = meshgrid(linspace(min(x_tr" << i + 1 << "),max(x_tr" << i + 1 << "),nx" << i + 1 << "),linspace(min(y_tr" << i + 1 << "),max(y_tr" << i + 1 << "),ny" << i + 1 << ")) ;\n";
            oss << "Z" << i + 1 << " =griddata(x_tr" << i + 1 << ",y_tr" << i + 1 << ",z_tr" << i + 1 << ",X" << i + 1 << ",Y" << i + 1 << ") ;\n";
            oss << "contour(X" << i + 1 << ",Y" << i + 1 << ",Z" << i + 1 << ") ;\n";
            oss << "hold on;\n";

            oss << "S_" << i + 1 << "= patch('Vertices', JXY_2Ds_" << i + 1 << ", 'Faces', JMs_" << i + 1 << ", 'FaceVertexCData', Data_" << i + 1 << ", 'FaceColor', 'interp', 'EdgeAlpha', 0.2, 'facealpha', 0);\n";
            oss << "hold on;\n";

            oss << "\n%%***flow rate vector 2d **************\n";
            oss << "q_vector_2d" << i + 1 << "=[";
            for (size_t j = 0; j < DFN_mesh.JXY[i].size(); ++j)
            {
                for (size_t yj = 0; yj < 3; ++yj)
                    oss << DFN_mesh.JXY[i][j](yj) << ", ";

                size_t Node_s = 0;
                if (DFN_mesh.Coe_Matr_guide[i][j](0) == -1)
                {
                    size_t i_frac = DFN_mesh.Coe_Matr_guide[i][j](1);
                    size_t j_node = DFN_mesh.Coe_Matr_guide[i][j](2);
                    Node_s = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                }
                else
                {
                    Node_s = DFN_mesh.Coe_Matr_guide[i][j](0);
                }

                oss << X_overall[Node_s] << ", ";
                oss << X_overall[Node_s + DFN_mesh.NO_all_pnts] << ", 0; ";
            }

            oss << "];\n";
            oss << "quiver3(q_vector_2d" << i + 1 << "(:,1), q_vector_2d" << i + 1 << "(:,2), q_vector_2d" << i + 1 << "(:,3), q_vector_2d" << i + 1 << "(:,4), q_vector_2d" << i + 1 << "(:,5), q_vector_2d" << i + 1 << "(:,6));\n\n";
            oss << "hold on;\n";

            oss << "xlabel('x (m)');\nylabel('y (m)');\ntitle('head contour and flow rate vector (Fig. " << i + 1 << ")');\nhold on;\n";
        }
        if (i == 0 && NO_frac != 1)
        {
            oss << "display('are you sure that you wanna show all fractures?');\n";
            oss << "pause;\n";
        }
    }
    oss.close();
};

inline void FEM_DFN::Verify_in_and_out_flux(DFN::DFN_mesh DFN_mesh, double *X_overall)
{
    double Q_in = 0, Q_out = 0;
    int NO_ele_in = 0, NO_ele_out = 0;
    for (size_t i = 0; i < DFN_mesh.Inlet.size(); ++i)
    {
        std::map<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>>::iterator inlet_A = DFN_mesh.Inlet[i].begin();

        while (inlet_A != DFN_mesh.Inlet[i].end())
        {
            size_t i_frac = inlet_A->first.first;
            size_t j_ele = inlet_A->first.second;
            size_t EdgNO = inlet_A->second.first(2);

            size_t pnt_f_0 = DFN_mesh.JM[i_frac][j_ele]((0 + EdgNO * 2) % 6);
            size_t pnt_f_1 = DFN_mesh.JM[i_frac][j_ele]((1 + EdgNO * 2) % 6);
            size_t pnt_f_2 = DFN_mesh.JM[i_frac][j_ele]((2 + EdgNO * 2) % 6);

            int NODE_f0 = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_0](0);
            int NODE_f1 = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_1](0);
            int NODE_f2 = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_2](0);

            if (NODE_f0 == -1)
            {
                size_t i_frac_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_0](1);
                size_t j_node_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_0](2);

                NODE_f0 = DFN_mesh.Coe_Matr_guide[i_frac_s][j_node_s](0);
            }

            if (NODE_f1 == -1)
            {
                size_t i_frac_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_1](1);
                size_t j_node_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_1](2);

                NODE_f1 = DFN_mesh.Coe_Matr_guide[i_frac_s][j_node_s](0);
            }

            if (NODE_f2 == -1)
            {
                size_t i_frac_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_2](1);
                size_t j_node_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_2](2);

                NODE_f2 = DFN_mesh.Coe_Matr_guide[i_frac_s][j_node_s](0);
            }

            double u_0 = X_overall[NODE_f0];
            double u_1 = X_overall[NODE_f1];
            double u_2 = X_overall[NODE_f2];

            double v_0 = X_overall[NODE_f0 + DFN_mesh.NO_all_pnts];
            double v_1 = X_overall[NODE_f1 + DFN_mesh.NO_all_pnts];
            double v_2 = X_overall[NODE_f2 + DFN_mesh.NO_all_pnts];

            double q0 = pow(pow(u_0, 2) + pow(v_0, 2), 0.5);
            double q1 = pow(pow(u_1, 2) + pow(v_1, 2), 0.5);
            double q2 = pow(pow(u_2, 2) + pow(v_2, 2), 0.5);
            Vector3d Vec_q_normal;
            Vec_q_normal << q0, q1, q2;
            cout << Vec_q_normal.transpose() << endl;

            double L = (DFN_mesh.JXY[i][pnt_f_0] - DFN_mesh.JXY[i][pnt_f_2]).norm();
            double Q_ele = Numerical_integration_linear_2(Vec_q_normal, L);
            Q_in += Q_ele;
            NO_ele_in++;
            inlet_A++;
        }
    };
    cout << "\n\n\n";
    for (size_t i = 0; i < DFN_mesh.Outlet.size(); ++i)
    {
        std::map<std::pair<size_t, size_t>, std::pair<Vector6d, Vector6d>>::iterator outlet_A = DFN_mesh.Outlet[i].begin();
        while (outlet_A != DFN_mesh.Outlet[i].end())
        {
            size_t i_frac = outlet_A->first.first;
            size_t j_ele = outlet_A->first.second;
            size_t EdgNO = outlet_A->second.first(2);

            size_t pnt_f_0 = DFN_mesh.JM[i_frac][j_ele]((0 + EdgNO * 2) % 6);
            size_t pnt_f_1 = DFN_mesh.JM[i_frac][j_ele]((1 + EdgNO * 2) % 6);
            size_t pnt_f_2 = DFN_mesh.JM[i_frac][j_ele]((2 + EdgNO * 2) % 6);

            int NODE_f0 = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_0](0);
            int NODE_f1 = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_1](0);
            int NODE_f2 = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_2](0);

            if (NODE_f0 == -1)
            {
                size_t i_frac_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_0](1);
                size_t j_node_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_0](2);

                NODE_f0 = DFN_mesh.Coe_Matr_guide[i_frac_s][j_node_s](0);
            }

            if (NODE_f1 == -1)
            {
                size_t i_frac_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_1](1);
                size_t j_node_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_1](2);

                NODE_f1 = DFN_mesh.Coe_Matr_guide[i_frac_s][j_node_s](0);
            }

            if (NODE_f2 == -1)
            {
                size_t i_frac_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_2](1);
                size_t j_node_s = DFN_mesh.Coe_Matr_guide[i_frac][pnt_f_2](2);

                NODE_f2 = DFN_mesh.Coe_Matr_guide[i_frac_s][j_node_s](0);
            }

            double u_0 = X_overall[NODE_f0];
            double u_1 = X_overall[NODE_f1];
            double u_2 = X_overall[NODE_f2];

            double v_0 = X_overall[NODE_f0 + DFN_mesh.NO_all_pnts];
            double v_1 = X_overall[NODE_f1 + DFN_mesh.NO_all_pnts];
            double v_2 = X_overall[NODE_f2 + DFN_mesh.NO_all_pnts];

            double q0 = pow(pow(u_0, 2) + pow(v_0, 2), 0.5);
            double q1 = pow(pow(u_1, 2) + pow(v_1, 2), 0.5);
            double q2 = pow(pow(u_2, 2) + pow(v_2, 2), 0.5);

            Vector3d Vec_q_normal;
            Vec_q_normal << q0, q1, q2;
            cout << Vec_q_normal.transpose() << endl;

            double L = (DFN_mesh.JXY[i][pnt_f_0] - DFN_mesh.JXY[i][pnt_f_2]).norm();
            double Q_ele = Numerical_integration_linear_2(Vec_q_normal, L);
            Q_out += Q_ele;
            NO_ele_out++;
            outlet_A++;
        }
    };

    cout << "Q_in: " << Q_in << ", "
         << "Q_out: " << Q_out << endl;
    cout << NO_ele_in << ", " << NO_ele_out << endl;
};

inline void FEM_DFN::Calculate_flux_of_A_Boun_ele(DFN::DFN_mesh DFN_mesh,
                                                  double *X_overall,
                                                  double K_coe,
                                                  size_t s_ifrac,
                                                  size_t s_jele,
                                                  size_t EdgeNO,
                                                  std::vector<std::pair<std::pair<Vector3d, Vector3d>, std::pair<Vector3d, Vector3d>>> &SF,
                                                  size_t NumOfPnt_q)
{
    Eigen::VectorXd h_e;
    h_e = Eigen::VectorXd::Zero(6);
    for (size_t yj = 0; yj < 6; ++yj)
    {
        size_t node_ID = DFN_mesh.JM[s_ifrac][s_jele](yj);
        int Id_guide = DFN_mesh.Coe_Matr_guide[s_ifrac][node_ID](0);
        if (Id_guide == -1)
        {
            size_t i_frac = DFN_mesh.Coe_Matr_guide[s_ifrac][node_ID](1);
            size_t j_node = DFN_mesh.Coe_Matr_guide[s_ifrac][node_ID](2);
            Id_guide = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
            if (Id_guide == -1)
            {
                cout << "Error! find wrong repetitive point!\n";
                exit(0);
            }
            h_e[yj] = X_overall[Id_guide];
        }
        else
        {
            h_e[yj] = X_overall[Id_guide];
        }
    }

    VectorXd xi_1D;
    xi_1D = Eigen::VectorXd::Zero(NumOfPnt_q);
    for (size_t ifr = 0; ifr < NumOfPnt_q; ++ifr)
    {
        xi_1D[ifr] = -1. + ifr * (2. / (NumOfPnt_q - 1.));
    };

    if (NumOfPnt_q != 3)
    {
        cout << "Error! Quadratic element!\n";
        exit(0);
    }

    for (size_t i = 0; i < NumOfPnt_q; ++i)
    {
        double xi, eta;
        //cout << "EdgeNO = " << EdgeNO << endl;
        if (EdgeNO == 0)
        {
            Vector3d RMK;
            Vector3d AS, BS;
            AS << 0, 0, 0;
            BS << 1, 0, 0;
            binary_linear_equation_Trans(AS, BS, RMK, xi_1D(i));

            xi = RMK(0);
            eta = RMK(1);
        }
        else if (EdgeNO == 1)
        {
            Vector3d RMK;
            Vector3d AS, BS;
            AS << 1, 0, 0;
            BS << 0, 1, 0;
            binary_linear_equation_Trans(AS, BS, RMK, xi_1D(i));

            xi = RMK(0);
            eta = RMK(1);
        }
        else if (EdgeNO == 2)
        {
            Vector3d RMK;
            Vector3d AS, BS;
            AS << 0, 1, 0;
            BS << 0, 0, 0;
            binary_linear_equation_Trans(AS, BS, RMK, xi_1D(i));

            xi = RMK(0);
            eta = RMK(1);
        }
        else
        {
            std::cout << "Edge NO of an element should be 0, 1 or 2 only!\n";
            exit(0);
        }
        xi = round(xi, 4);
        eta = round(eta, 4);
        //cout << " edge no " << EdgeNO << endl;
        //cout << " xi " << xi << " eta " << eta << endl;
        Vector3d pnt0 = DFN_mesh.JXY[s_ifrac][DFN_mesh.JM[s_ifrac][s_jele](0)], pnt2 = DFN_mesh.JXY[s_ifrac][DFN_mesh.JM[s_ifrac][s_jele](2)], pnt4 = DFN_mesh.JXY[s_ifrac][DFN_mesh.JM[s_ifrac][s_jele](4)];
        Vector3d pnt1 = DFN_mesh.JXY[s_ifrac][DFN_mesh.JM[s_ifrac][s_jele](1)], pnt3 = DFN_mesh.JXY[s_ifrac][DFN_mesh.JM[s_ifrac][s_jele](3)], pnt5 = DFN_mesh.JXY[s_ifrac][DFN_mesh.JM[s_ifrac][s_jele](5)];

        Eigen::VectorXd x_e, y_e;
        x_e = Eigen::VectorXd::Zero(6);
        y_e = x_e;
        x_e << pnt0(0), pnt1(0), pnt2(0), pnt3(0), pnt4(0), pnt5(0);
        y_e << pnt0(1), pnt1(1), pnt2(1), pnt3(1), pnt4(1), pnt5(1);
        double qx_2d;
        double qy_2d;

        Calculated_velocity_of_a_pnt_local(xi, eta, K_coe, h_e, x_e, y_e, qx_2d, qy_2d);

        double A_ele = Heron_formula(pnt0, pnt2, pnt4);

        if (i == 0 || i == 2)
        {
            Weighted_flux_shared_pnt(DFN_mesh, X_overall, K_coe, A_ele, qx_2d, qy_2d, s_ifrac, s_jele, EdgeNO, i);
        };
        //cout << "qx_final: " << qx_2d << ", qy_final: " << qy_2d << endl;
        //cout << endl;

        Vector3d Pnt_t_2d, q_2d, Pnt_t_3d, q_3d;
        Calculated_pnt_from_local_to_golbal(xi, eta, x_e, y_e, Pnt_t_2d);
        Pnt_t_3d = Pnt_t_2d;
        q_2d << qx_2d, qy_2d, 0;
        q_3d = q_2d;

        //-----rotation to 3D
        double R_angle_temp1 = DFN_mesh.Rota_angle[s_ifrac].first(0);
        Vector3d Frac_center;
        Frac_center << DFN_mesh.Rota_angle[s_ifrac].first(1), DFN_mesh.Rota_angle[s_ifrac].first(2), DFN_mesh.Rota_angle[s_ifrac].first(3);
        Vector3d temp3;
        temp3 << DFN_mesh.Rota_angle[s_ifrac].second(3), DFN_mesh.Rota_angle[s_ifrac].second(4), DFN_mesh.Rota_angle[s_ifrac].second(5);
        Vector3d Normal_frac;
        Normal_frac << DFN_mesh.Rota_angle[s_ifrac].second(0), DFN_mesh.Rota_angle[s_ifrac].second(1), DFN_mesh.Rota_angle[s_ifrac].second(2);

        if (abs(Normal_frac(0)) < 0.0001 && abs(Normal_frac(1)) < 0.0001 && abs(Normal_frac(2)) < 0.0001)
        {
            //nothing to do
            Pnt_t_3d += Frac_center;
        }
        else
        {
            Quaternion_t Q_axis_1;
            NormalizeRotation(R_angle_temp1, temp3, Q_axis_1);

            Vector3d temp4;
            Rotation(q_3d, Q_axis_1, temp4);
            q_3d = temp4;

            Vector3d temp5;
            Rotation(Pnt_t_3d, Q_axis_1, temp5);
            Pnt_t_3d = temp5 + Frac_center;
        }

        std::pair<std::pair<Vector3d, Vector3d>, std::pair<Vector3d, Vector3d>> RT;
        RT.first.first = Pnt_t_2d;
        RT.first.second = q_2d;
        RT.second.first = Pnt_t_3d;
        RT.second.second = q_3d;
        SF[i] = RT;
    }
};

inline void FEM_DFN::Weighted_flux_shared_pnt(DFN::DFN_mesh DFN_mesh,
                                              double *X_overall,
                                              double K_coe,
                                              double A_ele,
                                              double &qx,
                                              double &qy,
                                              size_t ty_ifrac,
                                              size_t ty_jele,
                                              size_t ty_edge_no,     // 0 - 2
                                              size_t ty_local_pnt_no // 0 - 2
)
{
    std::vector<Vector2d> overlap_q(1);
    overlap_q[0](0) = qx;
    overlap_q[0](1) = qy;
    std::vector<double> A_area(1);
    A_area[0] = A_ele;

    size_t golbal_pntID = DFN_mesh.JM[ty_ifrac][ty_jele]((ty_edge_no * (3 - 1) + ty_local_pnt_no) % 6);
    std::map<size_t, Vector2d> AAD;
    Vector2d RTY;
    RTY << ty_edge_no, ty_local_pnt_no;
    AAD.insert(std::make_pair(ty_jele, RTY));
    //ele no, edge no, local pnt no

    Recursion_find_ele_shared_a_pnt(DFN_mesh, ty_ifrac, ty_jele, golbal_pntID, ty_local_pnt_no, AAD);
    /*
    std::map<size_t, Vector2d>::iterator it_sy;
    it_sy = AAD.begin();
    while (it_sy != AAD.end())
    {
        cout << "neigh ele: " << it_sy->first << ", edge no: " << it_sy->second(0);
        cout << ", local pnt no: " << it_sy->second(1) << endl;
        it_sy++;
    }*/

    std::map<size_t, Vector2d>::iterator it_sy;
    it_sy = AAD.begin();
    while (it_sy != AAD.end())
    {
        //double qx_u = 0, qy_u = 0;
        size_t neigh_ele = it_sy->first;
        size_t edge_no_u = it_sy->second(0);
        size_t local_pnt_no = it_sy->second(1);
        //cout << "neigh ele: " << it_sy->first << ", edge no: " << it_sy->second(0);
        //cout << ", local pnt no: " << it_sy->second(1) << endl;
        if (neigh_ele == ty_jele)
        {
            //cout << "\tqx = " << overlap_q[0](0) << ", qy = " << overlap_q[0](1) << "; area: " << A_ele << "\n";
            it_sy++;
            continue;
        }

        Eigen::VectorXd h_e;
        h_e = Eigen::VectorXd::Zero(6);
        for (size_t yj = 0; yj < 6; ++yj)
        {
            size_t node_ID = DFN_mesh.JM[ty_ifrac][neigh_ele](yj);
            int Id_guide = DFN_mesh.Coe_Matr_guide[ty_ifrac][node_ID](0);
            if (Id_guide == -1)
            {
                size_t i_frac = DFN_mesh.Coe_Matr_guide[ty_ifrac][node_ID](1);
                size_t j_node = DFN_mesh.Coe_Matr_guide[ty_ifrac][node_ID](2);
                Id_guide = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
                if (Id_guide == -1)
                {
                    cout << "Error! find wrong repetitive point!\n";
                    exit(0);
                }
                h_e[yj] = X_overall[Id_guide];
            }
            else
            {
                h_e[yj] = X_overall[Id_guide];
            }
        }

        double xi, eta;

        if (edge_no_u == 0)
        {
            if (local_pnt_no == 0)
            {
                xi = 0;
                eta = 0;
            }
            else if (local_pnt_no == 2)
            {
                xi = 1;
                eta = 0;
            }
        }
        else if (edge_no_u == 1)
        {
            if (local_pnt_no == 0)
            {
                xi = 1;
                eta = 0;
            }
            else if (local_pnt_no == 2)
            {
                xi = 0;
                eta = 1;
            }
        }
        else if (edge_no_u == 2)
        {
            if (local_pnt_no == 0)
            {
                xi = 0;
                eta = 1;
            }
            else if (local_pnt_no == 2)
            {
                xi = 0;
                eta = 0;
            }
        }

        Vector3d pnt0 = DFN_mesh.JXY[ty_ifrac][DFN_mesh.JM[ty_ifrac][neigh_ele](0)], pnt2 = DFN_mesh.JXY[ty_ifrac][DFN_mesh.JM[ty_ifrac][neigh_ele](2)], pnt4 = DFN_mesh.JXY[ty_ifrac][DFN_mesh.JM[ty_ifrac][neigh_ele](4)];
        Vector3d pnt1 = DFN_mesh.JXY[ty_ifrac][DFN_mesh.JM[ty_ifrac][neigh_ele](1)], pnt3 = DFN_mesh.JXY[ty_ifrac][DFN_mesh.JM[ty_ifrac][neigh_ele](3)], pnt5 = DFN_mesh.JXY[ty_ifrac][DFN_mesh.JM[ty_ifrac][neigh_ele](5)];

        Eigen::VectorXd x_e, y_e;
        x_e = Eigen::VectorXd::Zero(6);
        y_e = x_e;
        x_e << pnt0(0), pnt1(0), pnt2(0), pnt3(0), pnt4(0), pnt5(0);
        y_e << pnt0(1), pnt1(1), pnt2(1), pnt3(1), pnt4(1), pnt5(1);
        double qx_2d = 0;
        double qy_2d = 0;
        Calculated_velocity_of_a_pnt_local(xi, eta, K_coe, h_e, x_e, y_e, qx_2d, qy_2d);

        double A_area_2 = Heron_formula(pnt0, pnt2, pnt4);

        Vector2d SUR;
        SUR << qx_2d, qy_2d;
        //cout << "\tqx = " << qx_2d << ", qy = " << qy_2d << "; area: " << A_area_2 << "\n";
        overlap_q.push_back(SUR);
        A_area.push_back(A_area_2);
        it_sy++;
    }

    qx = 0;
    qy = 0;

    double A_sum = 0;
    for (size_t i = 0; i < overlap_q.size(); ++i)
    {

        qx += overlap_q[i](0) * A_area[i];
        qy += overlap_q[i](1) * A_area[i];
        A_sum += A_area[i];
        //cout << qx << ", " << qy << ", " << A_sum << endl;
    }

    qx = qx / A_sum;
    qy = qy / A_sum;
};

inline void FEM_DFN::Recursion_find_ele_shared_a_pnt(DFN::DFN_mesh DFN_mesh,
                                                     size_t ifrac,
                                                     size_t jele,
                                                     size_t pnt_ID_p,
                                                     size_t local_pnt_ID_p,
                                                     std::map<size_t, Vector2d> &AAD)
{
    //we need a tag to show the computer that
    //if this element has been checked or not
    size_t IDX_yu = 0;
    for (size_t i = 0; i < DFN_mesh.neigh_shared[ifrac][jele].size(); ++i)
    {
        size_t neigh_ele = DFN_mesh.neigh_shared[ifrac][jele][i].first;
        size_t edge_no = DFN_mesh.neigh_shared[ifrac][jele][i].second(1);
        for (size_t j = 0; j < 3; j += 2)
        {
            size_t neigh_pnt = DFN_mesh.JM[ifrac][neigh_ele]((edge_no * 2 + j) % 6);
            if (neigh_pnt == pnt_ID_p)
            {
                Vector2d SG;
                SG << edge_no, j;
                std::pair<std::map<size_t, Vector2d>::iterator, bool> FS;
                FS = AAD.insert(std::make_pair(neigh_ele, SG));
                IDX_yu++;

                if (FS.second)
                    Recursion_find_ele_shared_a_pnt(DFN_mesh, ifrac, neigh_ele, pnt_ID_p,
                                                    j, AAD);
            }
        }
        if (i == DFN_mesh.neigh_shared[ifrac][jele].size() - 1 && IDX_yu == 0)
        {
            return;
        }
    }
};

inline void FEM_DFN::Calculated_velocity_of_a_pnt_local(double xi, double eta, double K_coe, VectorXd h_e, VectorXd x_e, VectorXd y_e, double &qx, double &qy)
{
    if (h_e.size() != 6 || x_e.size() != 6 || y_e.size() != 6)
    {
        cout << "in function 'Calculated_velocity_of_a_pnt_local', vector does not initialized!\n";
        exit(0);
    }

    Eigen::RowVectorXd p_N_p_x, p_N_p_y;
    p_N_p_x = Eigen::RowVectorXd::Zero(6);
    //p_N_p_x = Eigen::RowVectorXd::Zero(3);
    p_N_p_y = p_N_p_x;

    Eigen::VectorXd pd_N_over_pd_xi, pd_N_over_pd_eta;
    pd_N_over_pd_xi = Eigen::VectorXd::Zero(6);
    pd_N_over_pd_eta = pd_N_over_pd_xi;

    double N_natural_a_1,
        N_natural_a_2,
        N_natural_a_3,
        N_natural_a_4,
        N_natural_a_5,
        N_natural_a_6;

    double N_natural_b_1,
        N_natural_b_2,
        N_natural_b_3,
        N_natural_b_4,
        N_natural_b_5,
        N_natural_b_6;

    p_PHI_over_p_xi_and_eta(xi,
                            eta,
                            N_natural_a_1,
                            N_natural_a_2,
                            N_natural_a_3,
                            N_natural_a_4,
                            N_natural_a_5,
                            N_natural_a_6,
                            N_natural_b_1,
                            N_natural_b_2,
                            N_natural_b_3,
                            N_natural_b_4,
                            N_natural_b_5,
                            N_natural_b_6);

    //1 - 4 - 2 -5 -3 - 6
    pd_N_over_pd_xi << N_natural_a_1, N_natural_a_4, N_natural_a_2, N_natural_a_5, N_natural_a_3, N_natural_a_6;

    //1 - 4 - 2 -5 -3 - 6
    pd_N_over_pd_eta << N_natural_b_1, N_natural_b_4, N_natural_b_2, N_natural_b_5, N_natural_b_3, N_natural_b_6;

    Eigen::VectorXd pd_x_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_xi, pd_y_over_pd_eta;
    pd_x_over_pd_xi = Eigen::VectorXd::Zero(6);
    pd_x_over_pd_eta = pd_x_over_pd_xi;
    pd_y_over_pd_xi = pd_x_over_pd_xi;
    pd_y_over_pd_eta = pd_x_over_pd_xi;

    pd_x_over_pd_xi = pd_N_over_pd_xi.transpose() * x_e;
    pd_x_over_pd_eta = pd_N_over_pd_eta.transpose() * x_e;

    pd_y_over_pd_xi = pd_N_over_pd_xi.transpose() * y_e;
    pd_y_over_pd_eta = pd_N_over_pd_eta.transpose() * y_e;

    Eigen::MatrixXd Jacobi;
    Jacobi = Eigen::MatrixXd::Zero(2, 2);
    Jacobi << pd_x_over_pd_xi, pd_y_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_eta;

    Eigen::MatrixXd tem_o;
    tem_o = Eigen::MatrixXd::Zero(2, pd_N_over_pd_xi.rows());
    for (size_t op = 0; op < (size_t)tem_o.cols(); ++op)
    {
        tem_o(0, op) = pd_N_over_pd_xi(op);
        tem_o(1, op) = pd_N_over_pd_eta(op);
    }

    Eigen::MatrixXd AAA;
    AAA = Jacobi.inverse() * tem_o;

    for (size_t op = 0; op < 6; ++op)
    {
        p_N_p_x(op) = AAA(0, op);
        p_N_p_y(op) = AAA(1, op);
    }

    qx = -K_coe * p_N_p_x * h_e; // qx = -delta h / delta x * k,  or say, qx = - ix * k
    qy = -K_coe * p_N_p_y * h_e;
};

inline void FEM_DFN::Calculated_pnt_from_local_to_golbal(double xi, double eta, VectorXd x_e, VectorXd y_e, Vector3d &Center)
{
    double N1;
    double N2;
    double N3;
    double N4;
    double N5;
    double N6;
    PHI_shape_function(xi, eta, N1, N2, N3, N4, N5, N6);
    Eigen::VectorXd N;
    N = Eigen::VectorXd::Zero(6);
    N << N1, N4, N2, N5, N3, N6;

    Center(0) = x_e.dot(N);
    Center(1) = y_e.dot(N);
    Center(2) = 0;
};

inline void FEM_DFN::p_PHI_over_p_x_and_y(double xi,
                                          double eta,
                                          MatrixXd JXYe_x,
                                          MatrixXd JXYe_y,
                                          VectorXd &pd_N_over_pd_x,
                                          VectorXd &pd_N_over_pd_y,
                                          MatrixXd &Jacobi)
{
    if (JXYe_x.rows() != 6 || JXYe_x.cols() != 1 ||
        JXYe_y.rows() != 6 || JXYe_y.cols() != 1)
    {
        cout << "Error! in function 'p_PHI_over_p_x_and_y', JXYe_x does not initialize\n";
        exit(0);
    };
    Eigen::VectorXd pd_N_over_pd_xi, pd_N_over_pd_eta;
    pd_N_over_pd_xi = Eigen::VectorXd::Zero(6);
    pd_N_over_pd_eta = pd_N_over_pd_xi;

    double N_natural_a_1_s,
        N_natural_a_2_s,
        N_natural_a_3_s,
        N_natural_a_4_s,
        N_natural_a_5_s,
        N_natural_a_6_s;

    double N_natural_b_1_s,
        N_natural_b_2_s,
        N_natural_b_3_s,
        N_natural_b_4_s,
        N_natural_b_5_s,
        N_natural_b_6_s;

    p_PHI_over_p_xi_and_eta(xi,
                            eta,
                            N_natural_a_1_s,
                            N_natural_a_2_s,
                            N_natural_a_3_s,
                            N_natural_a_4_s,
                            N_natural_a_5_s,
                            N_natural_a_6_s,
                            N_natural_b_1_s,
                            N_natural_b_2_s,
                            N_natural_b_3_s,
                            N_natural_b_4_s,
                            N_natural_b_5_s,
                            N_natural_b_6_s);

    //1 - 4 - 2 -5 -3 - 6
    pd_N_over_pd_xi << N_natural_a_1_s, N_natural_a_4_s, N_natural_a_2_s, N_natural_a_5_s, N_natural_a_3_s, N_natural_a_6_s;

    pd_N_over_pd_eta << N_natural_b_1_s, N_natural_b_4_s, N_natural_b_2_s, N_natural_b_5_s, N_natural_b_3_s, N_natural_b_6_s;

    Eigen::VectorXd pd_x_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_xi, pd_y_over_pd_eta;
    pd_x_over_pd_xi = Eigen::VectorXd::Zero(6);
    pd_x_over_pd_eta = pd_x_over_pd_xi;
    pd_y_over_pd_xi = pd_x_over_pd_xi;
    pd_y_over_pd_eta = pd_x_over_pd_xi;

    pd_x_over_pd_xi = pd_N_over_pd_xi.transpose() * JXYe_x;
    pd_x_over_pd_eta = pd_N_over_pd_eta.transpose() * JXYe_x;

    pd_y_over_pd_xi = pd_N_over_pd_xi.transpose() * JXYe_y;
    pd_y_over_pd_eta = pd_N_over_pd_eta.transpose() * JXYe_y;

    Jacobi << pd_x_over_pd_xi, pd_y_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_eta;

    Eigen::MatrixXd tem_o;
    tem_o = Eigen::MatrixXd::Zero(2, pd_N_over_pd_xi.rows());
    for (size_t op = 0; op < (size_t)tem_o.cols(); ++op)
    {
        tem_o(0, op) = pd_N_over_pd_xi(op);
        tem_o(1, op) = pd_N_over_pd_eta(op);
    }

    Eigen::MatrixXd AAA;
    AAA = Jacobi.inverse() * tem_o;

    for (size_t op = 0; op < (size_t)pd_N_over_pd_x.rows(); ++op)
    {
        pd_N_over_pd_x(op) = AAA(0, op);
        pd_N_over_pd_y(op) = AAA(1, op);
    }
};

inline void FEM_DFN::PHI_shape_function(double xi,
                                        double eta,
                                        double &N_1,
                                        double &N_2,
                                        double &N_3,
                                        double &N_4,
                                        double &N_5,
                                        double &N_6)
{
    N_1 = (1 - xi - eta) * (1 - 2 * xi - 2 * eta);
    N_2 = xi * (2 * xi - 1);
    N_3 = eta * (2 * eta - 1);
    N_4 = 4 * xi * (1 - xi - eta);
    N_5 = 4 * xi * eta;
    N_6 = 4 * eta * (1 - xi - eta);
};

inline void FEM_DFN::p_PHI_over_p_xi_and_eta(double xi,
                                             double eta,
                                             double &N_natural_a_1_s,
                                             double &N_natural_a_2_s,
                                             double &N_natural_a_3_s,
                                             double &N_natural_a_4_s,
                                             double &N_natural_a_5_s,
                                             double &N_natural_a_6_s,
                                             double &N_natural_b_1_s,
                                             double &N_natural_b_2_s,
                                             double &N_natural_b_3_s,
                                             double &N_natural_b_4_s,
                                             double &N_natural_b_5_s,
                                             double &N_natural_b_6_s)
{
    N_natural_a_1_s = 4 * xi + 4 * eta - 3;
    N_natural_a_2_s = 4 * xi - 1;
    N_natural_a_3_s = 0;
    N_natural_a_4_s = 4 - 8 * xi - 4 * eta;
    N_natural_a_5_s = 4 * eta;
    N_natural_a_6_s = -4 * eta;

    N_natural_b_1_s = 4 * xi + 4 * eta - 3;
    N_natural_b_2_s = 0;
    N_natural_b_3_s = 4 * eta - 1;
    N_natural_b_4_s = -4 * xi;
    N_natural_b_5_s = 4 * xi;
    N_natural_b_6_s = 4 - 8 * eta - 4 * xi;
};

inline void FEM_DFN::p_PSI_over_p_x_and_y(double xi,
                                          double eta,
                                          MatrixXd JXYe_x,
                                          MatrixXd JXYe_y,
                                          VectorXd &pd_N_over_pd_x,
                                          VectorXd &pd_N_over_pd_y,
                                          MatrixXd &Jacobi)
{
    if (JXYe_x.rows() != 3 || JXYe_x.cols() != 1 ||
        JXYe_y.rows() != 3 || JXYe_y.cols() != 1)
    {
        cout << "Error! in function 'Calculated_pd_N_over_pd_x', JXYe_x does not initialize\n";
        exit(0);
    };
    Eigen::VectorXd pd_N_over_pd_xi, pd_N_over_pd_eta;
    pd_N_over_pd_xi = Eigen::VectorXd::Zero(3);
    pd_N_over_pd_eta = pd_N_over_pd_xi;

    double N_natural_a_1_s,
        N_natural_a_2_s,
        N_natural_a_3_s;

    double N_natural_b_1_s,
        N_natural_b_2_s,
        N_natural_b_3_s;

    p_PSI_over_p_xi_and_eta(xi,
                            eta,
                            N_natural_a_1_s,
                            N_natural_a_2_s,
                            N_natural_a_3_s,
                            N_natural_b_1_s,
                            N_natural_b_2_s,
                            N_natural_b_3_s);

    //1 - 2 - 3
    pd_N_over_pd_xi << N_natural_a_1_s, N_natural_a_2_s, N_natural_a_3_s;

    pd_N_over_pd_eta << N_natural_b_1_s, N_natural_b_2_s, N_natural_b_3_s;

    Eigen::VectorXd pd_x_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_xi, pd_y_over_pd_eta;
    pd_x_over_pd_xi = Eigen::VectorXd::Zero(3);
    pd_x_over_pd_eta = pd_x_over_pd_xi;
    pd_y_over_pd_xi = pd_x_over_pd_xi;
    pd_y_over_pd_eta = pd_x_over_pd_xi;

    pd_x_over_pd_xi = pd_N_over_pd_xi.transpose() * JXYe_x;
    pd_x_over_pd_eta = pd_N_over_pd_eta.transpose() * JXYe_x;

    pd_y_over_pd_xi = pd_N_over_pd_xi.transpose() * JXYe_y;
    pd_y_over_pd_eta = pd_N_over_pd_eta.transpose() * JXYe_y;

    Jacobi << pd_x_over_pd_xi, pd_y_over_pd_xi, pd_x_over_pd_eta, pd_y_over_pd_eta;

    Eigen::MatrixXd tem_o;
    tem_o = Eigen::MatrixXd::Zero(2, pd_N_over_pd_xi.rows());
    for (size_t op = 0; op < (size_t)tem_o.cols(); ++op)
    {
        tem_o(0, op) = pd_N_over_pd_xi(op);
        tem_o(1, op) = pd_N_over_pd_eta(op);
    }

    Eigen::MatrixXd AAA;
    AAA = Jacobi.inverse() * tem_o;

    for (size_t op = 0; op < (size_t)pd_N_over_pd_x.rows(); ++op)
    {
        pd_N_over_pd_x(op) = AAA(0, op);
        pd_N_over_pd_y(op) = AAA(1, op);
    }
};

inline void FEM_DFN::p_PSI_over_p_xi_and_eta(double xi,
                                             double eta,
                                             double &N_natural_a_1_s,
                                             double &N_natural_a_2_s,
                                             double &N_natural_a_3_s,
                                             double &N_natural_b_1_s,
                                             double &N_natural_b_2_s,
                                             double &N_natural_b_3_s)
{
    N_natural_a_1_s = -1;
    N_natural_a_2_s = 1;
    N_natural_a_3_s = 0;

    N_natural_b_1_s = -1;
    N_natural_b_2_s = 0;
    N_natural_b_3_s = 1;
};

inline void FEM_DFN::PSI_shape_function(double xi,
                                        double eta,
                                        double &N_1,
                                        double &N_2,
                                        double &N_3)
{
    N_1 = 1 - xi - eta;
    N_2 = xi;
    N_3 = eta;
};

inline void FEM_DFN::Three_D_velocity(DFN::DFN_mesh DFN_mesh, double *X_overall)
{
    Velo_3D.resize(DFN_mesh.JXY.size());
    for (size_t i = 0; i < DFN_mesh.JXY.size(); ++i)
    {
        Velo_3D[i].resize(DFN_mesh.JXY[i].size());
        for (size_t j = 0; j < DFN_mesh.JXY[i].size(); ++j)
        {
            int NODE_u = DFN_mesh.Coe_Matr_guide[i][j](0);
            if (NODE_u == -1)
            {
                size_t i_frac = DFN_mesh.Coe_Matr_guide[i][j](1);
                size_t j_node = DFN_mesh.Coe_Matr_guide[i][j](2);
                NODE_u = DFN_mesh.Coe_Matr_guide[i_frac][j_node](0);
            }

            double u = X_overall[NODE_u];
            double v = X_overall[NODE_u + DFN_mesh.NO_all_pnts];
            Velo_3D[i][j] << u, v, 0;

            //-----rotation to 3D
            double R_angle_temp1 = DFN_mesh.Rota_angle[i].first(0);
            Vector3d temp3;
            temp3 << DFN_mesh.Rota_angle[i].second(3), DFN_mesh.Rota_angle[i].second(4), DFN_mesh.Rota_angle[i].second(5);
            Vector3d Normal_frac;
            Normal_frac << DFN_mesh.Rota_angle[i].second(0), DFN_mesh.Rota_angle[i].second(1), DFN_mesh.Rota_angle[i].second(2);

            if (abs(Normal_frac(0)) < 0.0001 && abs(Normal_frac(1)) < 0.0001 && abs(Normal_frac(2)) < 0.0001)
            {
                //nothing
            }
            else
            {
                Quaternion_t Q_axis_1;
                NormalizeRotation(R_angle_temp1, temp3, Q_axis_1);

                Vector3d temp4;
                Rotation(Velo_3D[i][j], Q_axis_1, temp4);
                Velo_3D[i][j] = temp4;
            }
        }
    }
}
}; // namespace DFN