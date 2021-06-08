#pragma once

#include "../Mesh_H/Mesh_DFN.h"

namespace DFN
{

class FEM_DFN_A
{
public:
    std::vector<std::vector<double>> u_velocity;
    std::vector<std::vector<double>> v_velocity;
    std::vector<std::vector<Vector3d>> Velocity_3D;
    std::vector<std::vector<double>> h_head;

    std::vector<std::vector<Vector3d>> BC_top;
    std::vector<std::vector<Vector3d>> BC_bot;
    double Q_in;
    double Q_out;
    double Permeability;

public:
    FEM_DFN_A(DFN::Mesh_DFN DFN_mesh, DFN::Domain dom);
    void Assemble_overall_matrix(DFN::Mesh_DFN DFN_mesh, double *K_overall, double *F_overall, DFN::Domain dom);

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

    void matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN DFN_mesh, double *F_overall);

    void FEM_results(DFN::Mesh_DFN DFN_mesh, double *F_overall, DFN::Domain dom);

    void In_and_out_flux(DFN::Mesh_DFN DFN_mesh, DFN::Domain dom);
};

inline FEM_DFN_A::FEM_DFN_A(DFN::Mesh_DFN DFN_mesh, DFN::Domain dom)
{
    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    size_t Matrix_D = (NUM_NODES_velocity * 2 + NUM_NODES_p);

    //cout << NUM_NODES_velocity << ", " << NUM_NODES_p << endl;

    double *K_overall = new double[Matrix_D * Matrix_D];
    if (K_overall == NULL)
    {
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'K_overall'!\n");
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
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'F_overall'!\n");
    }
    else
    {
        for (size_t i = 0; i < Matrix_D; ++i)
        {
            F_overall[i] = 0;
        }
    }

    this->Assemble_overall_matrix(DFN_mesh, K_overall, F_overall, dom);
    // K x = F
    int n = Matrix_D; // dimensions of coefficient matrix K (2D)
    int m = 1;        // dimensions of F matrix
    //double a = 1, b = -1.00;

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
    }
    */

    //dgemv_("N", &n, &n, &a, K_overall, &n, F_overall, &m, &b, X_overall, &m);
    int *ipiv = new int[n];
    if (ipiv == NULL)
    {
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'F_overall'!\n");
    }
    else
    {
        for (size_t i = 0; i < (size_t)n; ++i)
        {
            ipiv[i] = 0;
        }
    }
    int info;
    std::cout << "\033[32mstart solving matrix;\n\033[0m";
    std::cout << "The size of matrix is " << Matrix_D << endl;
    dgesv_(&n, &m, K_overall, &n, ipiv, F_overall, &n, &info);
    std::cout << "\033[32mfinish solving matrix;\n\033[0m";;
    // now F_overall is the solution X
    if (info > 0)
    {
        throw Error_throw_ignore("The solution could not be computed!!\n");
    }
    /*
    cout << "\nX_overall with BLAS;\n";
    for (size_t i = 0; i < NUM_NODES_velocity; ++i)
    {
        cout << F_overall[i] << endl;
    }
    */

    this->FEM_results(DFN_mesh, F_overall, dom);

    //this->matlab_plot("FEM_DFN.mat", "FEM_DFN.m", dom, DFN_mesh, F_overall);

    this->In_and_out_flux(DFN_mesh, dom);
    //cout << "in : " << this->Q_in << endl;
    //cout << "out: " << this->Q_out << endl;

    delete[] ipiv;
    ipiv = NULL;
    delete[] F_overall;
    F_overall = NULL;
    delete[] K_overall;
    K_overall = NULL;
};

inline void FEM_DFN_A::Assemble_overall_matrix(DFN::Mesh_DFN DFN_mesh, double *K_overall, double *F_overall, DFN::Domain dom)
{
    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    size_t Matrix_D = (NUM_NODES_velocity * 2 + NUM_NODES_p);

    for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
    {
        size_t Frac_Tag = DFN_mesh.Frac_Tag[i];

        double Kper = dom.Fractures[Frac_Tag].Conductivity;

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

                    this->PHI_shape_function(xi(ik),
                                             eta(jk),
                                             Phi_1,
                                             Phi_2,
                                             Phi_3,
                                             Phi_4,
                                             Phi_5,
                                             Phi_6);

                    PHI << Phi_1, Phi_4, Phi_2, Phi_5, Phi_3, Phi_6;

                    this->p_PHI_over_p_x_and_y(xi(ik), eta(jk), JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

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
                    this->PHI_shape_function(xi(ik),
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

                    this->p_PSI_over_p_x_and_y(xi(ik), eta(jk), JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

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

            //------------------ assemble
            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    size_t Node_m = DFN_mesh.JM[i][j](jq);
                    size_t Node_n = DFN_mesh.JM[i][j](kq);

                    if (DFN_mesh.Guide_frac[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Guide_frac[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](kq)](0);
                    }

                    size_t Global_Idx = Node_m * Matrix_D + Node_n;

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(D1e(jq, kq), 4);
                }
            }

            //
            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    size_t Node_m = DFN_mesh.JM[i][j](jq);
                    size_t Node_n = DFN_mesh.JM[i][j](kq);
                    if (DFN_mesh.Guide_frac[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Guide_frac[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](kq)](0);
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
                    if (DFN_mesh.Guide_frac[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Guide_frac[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Guide_frac[i_frac][j_node](3);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](yt)](3);
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
                    if (DFN_mesh.Guide_frac[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](jq)](0);
                    }

                    if (DFN_mesh.Guide_frac[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Guide_frac[i_frac][j_node](3);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](yt)](3);
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
                    if (DFN_mesh.Guide_frac[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Guide_frac[i_frac][j_node](3);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](yt)](3);
                    }

                    if (DFN_mesh.Guide_frac[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](kq)](0);
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
                    if (DFN_mesh.Guide_frac[i][Node_m](0) == -1) // repetitive point
                    {
                        //cout << "find -1;\n";
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_m](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_m](2); // the j th node
                        //cout << "i_frac " << i_frac << ", j_node: " << j_node << endl;
                        Node_m = DFN_mesh.Guide_frac[i_frac][j_node](3);
                    }
                    else
                    {
                        Node_m = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](yt)](3);
                    }

                    if (DFN_mesh.Guide_frac[i][Node_n](0) == -1) // repetitive point
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][Node_n](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][Node_n](2); // the j th node
                        Node_n = DFN_mesh.Guide_frac[i_frac][j_node](0);
                    }
                    else
                    {
                        Node_n = DFN_mesh.Guide_frac[i][DFN_mesh.JM[i][j](kq)](0);
                    }

                    size_t Global_Idx = (Node_m + 2 * NUM_NODES_velocity) * Matrix_D + (Node_n + NUM_NODES_velocity);

                    K_overall[Global_Idx] = K_overall[Global_Idx] + round(B2e(jq, kq), 4);
                }
            }
        }
    }

    // boundary condition
    // i = frac
    // j = element

    this->BC_top.resize(DFN_mesh.JM.size());
    this->BC_bot.resize(DFN_mesh.JM.size());
    //2nd bc

    for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
    {
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            for (size_t ik = 0; ik < 6; ik += 2)
            {
                size_t edge_0 = ik;
                size_t edge_1 = ik + 1;
                size_t edge_2 = (ik + 2) % 6;

                size_t pnt_ID_0 = DFN_mesh.JM[i][j](edge_0);
                size_t pnt_ID_1 = DFN_mesh.JM[i][j](edge_1);
                size_t pnt_ID_2 = DFN_mesh.JM[i][j](edge_2);

                //-----------find if there are top or bottom edge-------
                if (DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_top == true &&
                    DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_top == true &&
                    DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_top == true)
                {
                    BC_top[i].push_back(Vector3d{(double)pnt_ID_0, (double)pnt_ID_1, (double)pnt_ID_2}); //
                }
                else if (DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_bottom == true &&
                         DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_bottom == true &&
                         DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_bottom == true)
                {
                    BC_bot[i].push_back(Vector3d{(double)pnt_ID_0, (double)pnt_ID_1, (double)pnt_ID_2}); //
                }

                // find lateral boundary
                if (((DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_top == false || DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_top == false || DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_top == false) &&
                     (DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_bottom == false || DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_bottom == false || DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_bottom == false)) &&
                    ((DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_front == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_front == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_front == true) ||
                     (DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_back == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_back == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_back == true) ||
                     (DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_left == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_left == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_left == true) ||
                     (DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_0].If_model_right == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_1].If_model_right == true && DFN_mesh.Pnt_attri_Mo_B[i][pnt_ID_2].If_model_right == true) ||
                     (DFN_mesh.Pnt_attribute[i][pnt_ID_0].If_frac_bound == true && DFN_mesh.Pnt_attribute[i][pnt_ID_1].If_frac_bound == true && DFN_mesh.Pnt_attribute[i][pnt_ID_2].If_frac_bound == true)))
                {

                    double q0 = 0, q1 = 0, q2 = 0; // velocity boundary condition

                    double L = (DFN_mesh.JXY[i][pnt_ID_0] - DFN_mesh.JXY[i][pnt_ID_2]).norm();
                    Vector2d Ge;
                    Ge(0) = (L * (q0 + 2 * q1)) / 3;
                    Ge(1) = (L * (2 * q1 + q2)) / 3;

                    size_t Node_s1 = 0, Node_s2 = 0;

                    if (DFN_mesh.Guide_frac[i][pnt_ID_0](0) == -1)
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][pnt_ID_0](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][pnt_ID_0](2); // the j th node
                        Node_s1 = DFN_mesh.Guide_frac[i_frac][j_node](3);
                    }
                    else
                    {
                        Node_s1 = DFN_mesh.Guide_frac[i][pnt_ID_0](3);
                    }

                    if (DFN_mesh.Guide_frac[i][pnt_ID_2](0) == -1)
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][pnt_ID_2](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][pnt_ID_2](2); // the j th node
                        Node_s2 = DFN_mesh.Guide_frac[i_frac][j_node](3);
                    }
                    else
                    {
                        Node_s2 = DFN_mesh.Guide_frac[i][pnt_ID_2](3);
                    }

                    F_overall[Node_s1 + NUM_NODES_velocity * 2] = F_overall[Node_s1 + NUM_NODES_velocity * 2] + Ge(0);

                    F_overall[Node_s2 + NUM_NODES_velocity * 2] = F_overall[Node_s2 + NUM_NODES_velocity * 2] + Ge(1);
                }
            }
        }
    }

    // 1st BC
    // i th fracture
    for (size_t i = 0; i < DFN_mesh.JXY.size(); ++i)
    {
        for (size_t j = 0; j < DFN_mesh.JXY[i].size(); ++j)
        {
            if (DFN_mesh.Guide_frac[i][j](3) != -2)
            {
                if (DFN_mesh.Pnt_attri_Mo_B[i][j].If_model_top == true || DFN_mesh.Pnt_attri_Mo_B[i][j].If_model_bottom == true)
                {
                    size_t node_s = 0;
                    if (DFN_mesh.Guide_frac[i][j](0) == -1)
                    {
                        size_t i_frac = DFN_mesh.Guide_frac[i][j](1); // the i th frac
                        size_t j_node = DFN_mesh.Guide_frac[i][j](2); // the j th node
                        node_s = DFN_mesh.Guide_frac[i_frac][j_node](3);
                    }
                    else
                    {
                        node_s = DFN_mesh.Guide_frac[i][j](3);
                    }
                    double BC_1 = 0;

                    if (DFN_mesh.Pnt_attri_Mo_B[i][j].If_model_top == true)
                    {

                        BC_1 = 100;
                    }

                    if (DFN_mesh.Pnt_attri_Mo_B[i][j].If_model_bottom == true)
                    {

                        BC_1 = 20;
                    }

                    for (size_t go = 0; go < Matrix_D; ++go)
                    {
                        size_t golbalID = go * Matrix_D + node_s + 2 * NUM_NODES_velocity;

                        F_overall[go] = F_overall[go] - K_overall[golbalID] * BC_1;

                        //cout << "K_overall[golbalID]: " << K_overall[golbalID] << endl;

                        K_overall[golbalID] = 0;

                        size_t golbalID_row = (node_s + 2 * NUM_NODES_velocity) * Matrix_D + go;
                        K_overall[golbalID_row] = 0;
                    }
                    //exit(0);
                    K_overall[(node_s + 2 * NUM_NODES_velocity) * Matrix_D + (node_s + 2 * NUM_NODES_velocity)] = 1;

                    F_overall[(node_s + 2 * NUM_NODES_velocity)] = BC_1;
                }
            }
        }
    }
}

inline void FEM_DFN_A::p_PHI_over_p_x_and_y(double xi,
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
        throw Error_throw_pause("Error! in function 'p_PHI_over_p_x_and_y', JXYe_x does not initialize\n");
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

inline void FEM_DFN_A::PHI_shape_function(double xi,
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

inline void FEM_DFN_A::p_PHI_over_p_xi_and_eta(double xi,
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

inline void FEM_DFN_A::p_PSI_over_p_x_and_y(double xi,
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
        throw Error_throw_pause("Error! in function 'Calculated_pd_N_over_pd_x', JXYe_x does not initialize\n");
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

inline void FEM_DFN_A::p_PSI_over_p_xi_and_eta(double xi,
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

inline void FEM_DFN_A::PSI_shape_function(double xi,
                                          double eta,
                                          double &N_1,
                                          double &N_2,
                                          double &N_3)
{
    N_1 = 1 - xi - eta;
    N_2 = xi;
    N_3 = eta;
};

inline void FEM_DFN_A::FEM_results(DFN::Mesh_DFN DFN_mesh, double *F_overall, DFN::Domain dom)
{
    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    //size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    //
    u_velocity.resize(DFN_mesh.JXY.size());
    v_velocity.resize(DFN_mesh.JXY.size());
    h_head.resize(DFN_mesh.JXY.size());
    Velocity_3D.resize(DFN_mesh.JXY.size());

    for (size_t i = 0; i < DFN_mesh.JXY.size(); ++i)
    {
        u_velocity[i].resize(DFN_mesh.JXY[i].size());
        v_velocity[i].resize(DFN_mesh.JXY[i].size());
        h_head[i].resize(DFN_mesh.JXY[i].size());
        Velocity_3D[i].resize(DFN_mesh.JXY[i].size());

        DFN::Polygon_convex_3D This_frac{dom.Fractures[DFN_mesh.Frac_Tag[i]].Verts};
        Vector3d temp1;
        DFN::Vector_2 v(This_frac.Normal_vector, temp1);
        double R_angle_temp1 = 0;
        double x_temp = This_frac.Beta;
        R_angle_temp1 = x_temp * M_PI / 180;
        Quaternion_t Q_axis_1;

        for (size_t j = 0; j < DFN_mesh.JXY[i].size(); ++j)
        {

            size_t pntTag = 0;
            if (DFN_mesh.Guide_frac[i][j](0) == -1)
            {
                size_t i_frac = DFN_mesh.Guide_frac[i][j](1); // the i th frac
                size_t j_node = DFN_mesh.Guide_frac[i][j](2); // the j th node
                pntTag = DFN_mesh.Guide_frac[i_frac][j_node](0);
            }
            else
            {
                pntTag = DFN_mesh.Guide_frac[i][j](0);
            }

            u_velocity[i][j] = F_overall[pntTag];
            v_velocity[i][j] = F_overall[pntTag + NUM_NODES_velocity];

            //-----------calculate 3D velocity-----------------------------------
            std::vector<Vector3d> Verts_1(1), Verts_2(1);
            Verts_1[0] << u_velocity[i][j], v_velocity[i][j], 0;

            if (x_temp > 0.0001)
            {
                DFN::Rotation_verts R1(Verts_1, -R_angle_temp1, Q_axis_1, temp1, Verts_2);
            }
            else
            {
                Verts_2[0] = Verts_1[0];
            }
            //--------------------------------------------------------------------
            Velocity_3D[i][j] = Verts_2[0];
            //cout << Verts_2[0].transpose() << endl;

            if (DFN_mesh.Guide_frac[i][j](3) != -2) // it is a corner point
            {
                pntTag = 0;
                if (DFN_mesh.Guide_frac[i][j](0) == -1)
                {
                    size_t i_frac = DFN_mesh.Guide_frac[i][j](1); // the i th frac
                    size_t j_node = DFN_mesh.Guide_frac[i][j](2); // the j th node
                    pntTag = DFN_mesh.Guide_frac[i_frac][j_node](3);
                }
                else
                {
                    pntTag = DFN_mesh.Guide_frac[i][j](3);
                }

                h_head[i][j] = F_overall[pntTag + 2 * NUM_NODES_velocity];
            }
            else
            {
                h_head[i][j] = 0;
            }
        }
    }
};

inline void FEM_DFN_A::In_and_out_flux(DFN::Mesh_DFN DFN_mesh, DFN::Domain dom)
{
    this->Q_in = 0;
    this->Q_out = 0;
    for (size_t i = 0; i < this->BC_top.size(); ++i)
    {
        for (size_t j = 0; j < this->BC_top[i].size(); ++j)
        {
            size_t pointID0 = this->BC_top[i][j](0);
            //size_t pointID1 = this->BC_top[i][j](1);
            size_t pointID2 = this->BC_top[i][j](2);

            double L = (DFN_mesh.JXY[i][pointID0] - DFN_mesh.JXY[i][pointID2]).norm();
            Vector2d V0 = Vector2d{this->u_velocity[i][pointID0], this->v_velocity[i][pointID0]};
            //Vector2d V1 = Vector2d{this->u_velocity[i][pointID1], this->v_velocity[i][pointID1]};
            Vector2d V2 = Vector2d{this->u_velocity[i][pointID2], this->v_velocity[i][pointID2]};
            double ave_velocity = (V0.norm() + V2.norm()) / 2.0;

            double K = 1.0;

            Q_in += L * K * ave_velocity;
        }
    }

    for (size_t i = 0; i < this->BC_bot.size(); ++i)
    {
        for (size_t j = 0; j < this->BC_bot[i].size(); ++j)
        {
            size_t pointID0 = this->BC_bot[i][j](0);
            //size_t pointID1 = this->BC_bot[i][j](1);
            size_t pointID2 = this->BC_bot[i][j](2);

            double L = (DFN_mesh.JXY[i][pointID0] - DFN_mesh.JXY[i][pointID2]).norm();
            Vector2d V0 = Vector2d{this->u_velocity[i][pointID0], this->v_velocity[i][pointID0]};
            //Vector2d V1 = Vector2d{this->u_velocity[i][pointID1], this->v_velocity[i][pointID1]};
            Vector2d V2 = Vector2d{this->u_velocity[i][pointID2], this->v_velocity[i][pointID2]};
            double ave_velocity = (V0.norm() + V2.norm()) / 2.0;

            double K = 1.0;

            Q_out += L * K * ave_velocity;
        }
    }

    if(abs(Q_in - Q_out) / (Q_in > Q_out? Q_in : Q_out) < 0.3)
    {
            double bc_top = 100, bc_bot = 20;
            double L = dom.Model_domain(0) - dom.Model_domain(1);

            double gradient_head = (bc_top - bc_bot) / L;
            //cout << "gradient_head: " << gradient_head << endl;
            //cout << "Q: " << (Q_out + Q_in) * 0.5 << endl;
            this->Permeability = (Q_out + Q_in) * 0.5 / (pow(L, 2.0) * gradient_head);
    }
    else
    {
        cout << "Q_out: " << Q_out << endl;
        cout << "Q_in: " << Q_in << endl;
        throw Error_throw_ignore("Found large difference between inlet and outlet flux!\n");
    }
};

inline void FEM_DFN_A::matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN DFN_mesh, double *F_overall)
{
    // size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    //size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;

    // mat file

    const char *filename = FileKey_mat.c_str();
    MATFile *pMatFile;
    pMatFile = matOpen(filename, "w");

    if (!pMatFile)
    {
        throw Error_throw_ignore("cannot create mat file in class Mesh_DFN\n");
    }
    // frac

    for (size_t i = 0; i < DFN_mesh.Frac_Tag.size(); ++i)
    {
        //cout << "i: " << i << endl;
        size_t len = dom.Fractures[DFN_mesh.Frac_Tag[i]].Verts_trim.size(); // number of verts

        double *pData1;
        double *pData2;
        double *pData3;
        double *pData4;
        double *pData5;
        double *pData6;
        double *pData7;
        double *pData8;
        double *pData9;

        pData1 = (double *)mxCalloc(len, sizeof(double));
        pData2 = (double *)mxCalloc(len, sizeof(double));
        pData3 = (double *)mxCalloc(len, sizeof(double));
        pData4 = (double *)mxCalloc(DFN_mesh.JXY_3D[i].size() * 3, sizeof(double));
        pData5 = (double *)mxCalloc(DFN_mesh.JM[i].size() * 6, sizeof(double));
        pData6 = (double *)mxCalloc(DFN_mesh.JXY[i].size(), sizeof(double));
        pData7 = (double *)mxCalloc(DFN_mesh.JXY[i].size(), sizeof(double));
        pData8 = (double *)mxCalloc(DFN_mesh.JXY[i].size(), sizeof(double));
        pData9 = (double *)mxCalloc(DFN_mesh.JXY[i].size() * 3, sizeof(double));

        mxArray *pMxArray1;
        mxArray *pMxArray2;
        mxArray *pMxArray3;
        mxArray *pMxArray4;
        mxArray *pMxArray5;
        mxArray *pMxArray6;
        mxArray *pMxArray7;
        mxArray *pMxArray8;
        mxArray *pMxArray9;

        pMxArray1 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray2 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray3 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray4 = mxCreateDoubleMatrix(DFN_mesh.JXY_3D[i].size(), 3, mxREAL);
        pMxArray5 = mxCreateDoubleMatrix(DFN_mesh.JM[i].size(), 6, mxREAL);
        pMxArray6 = mxCreateDoubleMatrix(DFN_mesh.JXY[i].size(), 1, mxREAL);
        pMxArray7 = mxCreateDoubleMatrix(DFN_mesh.JXY[i].size(), 1, mxREAL);
        pMxArray8 = mxCreateDoubleMatrix(DFN_mesh.JXY[i].size(), 1, mxREAL);
        pMxArray9 = mxCreateDoubleMatrix(DFN_mesh.JXY[i].size(), 3, mxREAL);

        if (!pMxArray1 || !pMxArray2 || !pMxArray3 || !pMxArray4 || !pMxArray5 || !pMxArray6 || !pMxArray7 || !pMxArray8 || !pMxArray9)
        {
            throw Error_throw_pause("cannot create pMxArray in class FEM_DFN_A\n");
        }

        if (!pData1 || !pData2 || !pData3 || !pData4 || !pData5 || !pData6 || !pData7 || !pData8 || !pData9)
        {
            throw Error_throw_pause("cannot create pData in class FEM_DFN_A\n");
        }

        for (size_t j = 0; j < len; j++)
        {
            pData1[j] = dom.Fractures[DFN_mesh.Frac_Tag[i]].Verts_trim[j](0);
            pData2[j] = dom.Fractures[DFN_mesh.Frac_Tag[i]].Verts_trim[j](1);
            pData3[j] = dom.Fractures[DFN_mesh.Frac_Tag[i]].Verts_trim[j](2);
        }

        for (size_t j = 0; j < DFN_mesh.JXY_3D[i].size() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / DFN_mesh.JXY_3D[i].size()); // column
            l = j % DFN_mesh.JXY_3D[i].size();       // row

            pData4[j] = DFN_mesh.JXY_3D[i][l](k);
        }

        for (size_t j = 0; j < DFN_mesh.JM[i].size() * 6; ++j)
        {
            size_t k, l;
            k = ceil(j / DFN_mesh.JM[i].size()); // column
            l = j % DFN_mesh.JM[i].size();       // row

            pData5[j] = DFN_mesh.JM[i][l](k) + 1;
        }

        for (size_t j = 0; j < DFN_mesh.JXY[i].size(); ++j)
        {
            pData6[j] = u_velocity[i][j];
            pData7[j] = v_velocity[i][j];

            pData8[j] = h_head[i][j];
        }

        for (size_t j = 0; j < DFN_mesh.JXY[i].size() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / DFN_mesh.JXY_3D[i].size()); // column
            l = j % DFN_mesh.JXY_3D[i].size();       // row
            pData9[j] = this->Velocity_3D[i][l](k);  // wrong!!!!!
        }

        mxSetData(pMxArray1, pData1);
        mxSetData(pMxArray2, pData2);
        mxSetData(pMxArray3, pData3);
        mxSetData(pMxArray4, pData4);
        mxSetData(pMxArray5, pData5);
        mxSetData(pMxArray6, pData6);
        mxSetData(pMxArray7, pData7);
        mxSetData(pMxArray8, pData8);
        mxSetData(pMxArray9, pData9);

        string ft = to_string(i);

        string Fracx = "Frac_" + ft + "_x";
        string Fracy = "Frac_" + ft + "_y";
        string Fracz = "Frac_" + ft + "_z";
        string FracJXY3D = "Frac_" + ft + "_JXY3D";
        string FracJM = "Frac_" + ft + "_JM";
        string FracU = "Frac_" + ft + "_u";
        string FracV = "Frac_" + ft + "_v";
        string FracH = "Frac_" + ft + "_h";
        string Frac3DV = "Frac_" + ft + "_3Dvelocity";

        const char *Fracx_s = Fracx.c_str();
        const char *Fracy_s = Fracy.c_str();
        const char *Fracz_s = Fracz.c_str();
        const char *FracJXY3D_s = FracJXY3D.c_str();
        const char *FracJM_s = FracJM.c_str();
        const char *FracU_s = FracU.c_str();
        const char *FracV_s = FracV.c_str();
        const char *FracH_s = FracH.c_str();
        const char *Frac3DV_s = Frac3DV.c_str();

        matPutVariable(pMatFile, Fracx_s, pMxArray1);
        matPutVariable(pMatFile, Fracy_s, pMxArray2);
        matPutVariable(pMatFile, Fracz_s, pMxArray3);
        matPutVariable(pMatFile, FracJXY3D_s, pMxArray4);
        matPutVariable(pMatFile, FracJM_s, pMxArray5);
        matPutVariable(pMatFile, FracU_s, pMxArray6);
        matPutVariable(pMatFile, FracV_s, pMxArray7);
        matPutVariable(pMatFile, FracH_s, pMxArray8);
        matPutVariable(pMatFile, Frac3DV_s, pMxArray9);

        mxFree(pData1);
        mxFree(pData2);
        mxFree(pData3);
        mxFree(pData4);
        mxFree(pData5);
        mxFree(pData6);
        mxFree(pData7);
        mxFree(pData8);
        mxFree(pData9);
    }

    matClose(pMatFile);

    // m file
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    oss << "%% matrix dimension is " << NUM_NODES_velocity * 2 + NUM_NODES_p << endl;
    oss << "load('" << FileKey_mat << "');\n";
    oss << "figure(1);\n";
    oss << "title('DFN');\n";
    for (size_t i = 0; i < DFN_mesh.Frac_Tag.size(); ++i)
    {
        oss << "fill3([Frac_" << i << "_x; Frac_" << i << "_x(1,1)], [Frac_" << i << "_y; Frac_" << i << "_y(1,1)], [Frac_" << i << "_z; Frac_" << i << "_z(1,1)], [rand rand rand]);\nhold on;\n";
    }
    oss << "\n\n";
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

    oss << "\n\nfigure(2);\n";
    oss << "view(3);\n";
    oss << "title('u velocity');\n";
    for (size_t i = 0; i < DFN_mesh.Frac_Tag.size(); ++i)
    {
        oss << "P" << i << " = patch('Vertices', Frac_" << i << "_JXY3D, 'Faces', Frac_" << i << "_JM, 'FaceVertexCData', Frac_" << i << "_u, 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1);\n";
        oss << "hold on;\n";
    }
    oss << "colorbar;\n";
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

    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\n";

    oss << "\n\nfigure(3);\n";
    oss << "view(3);\n";
    oss << "title('v velocity');\n";
    for (size_t i = 0; i < DFN_mesh.Frac_Tag.size(); ++i)
    {
        oss << "P" << i << " = patch('Vertices', Frac_" << i << "_JXY3D, 'Faces', Frac_" << i << "_JM, 'FaceVertexCData', Frac_" << i << "_v, 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1);\n";
        oss << "hold on;\n";
    }
    oss << "colorbar;\n";
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

    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\n";

    oss << "\n\nfigure(4);\n";
    oss << "view(3);\n";
    oss << "title('head and 3D velocity vector');\n";
    for (size_t i = 0; i < DFN_mesh.Frac_Tag.size(); ++i)
    {
        oss << "P" << i << " = patch('Vertices', Frac_" << i << "_JXY3D, 'Faces', Frac_" << i << "_JM(:, [1,3,5]), 'FaceVertexCData', Frac_" << i << "_h, 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1);\n";
        oss << "hold on;\n";
    }

    oss << "all_points = [\n";
    for (size_t i = 0; i < DFN_mesh.Frac_Tag.size(); ++i)
    {
        oss << "\tFrac_" << i << "_JXY3D;\n";
    }
    oss << "];\n";

    oss << "all_3D_velocity = [\n";
    for (size_t i = 0; i < DFN_mesh.Frac_Tag.size(); ++i)
    {
        oss << "Frac_" << i << "_3Dvelocity;\n";
    }
    oss << "];\n";
    oss << "quiver3(all_points(:,1), all_points(:,2), all_points(:,3), all_3D_velocity(:,1), all_3D_velocity(:,2), all_3D_velocity(:,3), 'r', 'linewidth', 1.2);\n\n";
    oss << "hold on;\n";
    oss << "colorbar;\n";
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
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\n";

    oss.close();
};

}; // namespace DFN
