#pragma once

//#include "../Mesh_H/Mesh_DFN.h"
#include "../Geometry_H/Normal_vector_2D.h"
#include "../Mesh_H/Mesh_DFN_overall.h"
#include "../Using_UMFPACK/Using_UMFPACK.h"
namespace DFN
{

class FEM_DFN_coupling
{
public:
    std::vector<Vector3d> Velocity_3D;
    VectorXd h_head;

    std::vector<std::vector<Vector3d>> BC_top;
    std::vector<std::vector<Vector3d>> BC_bot;
    double Q_in;
    double Q_out;
    double Permeability = 0;

    double *F_overall;

public:
    FEM_DFN_coupling(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom);
    void Assemble_overall_matrix(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, double *F_overall, DFN::Domain dom);

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

    void matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN_overall DFN_mesh, double *F_overall);

    void FEM_results(DFN::Mesh_DFN_overall DFN_mesh, double *F_overall, DFN::Domain dom);

    void In_and_out_flux(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom);

    ~FEM_DFN_coupling();
};

inline FEM_DFN_coupling::FEM_DFN_coupling(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom)
{

    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    size_t Matrix_D = (NUM_NODES_velocity * 2 + NUM_NODES_p);

    if (Matrix_D > 4e4)
    {
        throw Error_throw_pause("The dimension of the matrix is too large!\n");
    }

    double *K_overall = new double[Matrix_D * Matrix_D];

    if (K_overall == NULL)
    {
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'K_overall'!\n");
    }

    F_overall = new double[Matrix_D];
    if (F_overall == NULL)
    {
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'F_overall'!\n");
    }

    this->Assemble_overall_matrix(DFN_mesh, K_overall, F_overall, dom);

    // K x = F
    //int n = Matrix_D; // dimensions of coefficient matrix K (2D)
    //int m = 1;        // dimensions of F matrix
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
    /*
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
    //std::cout << "The size of matrix is " << Matrix_D << endl;
    dgesv_(&n, &m, K_overall, &n, ipiv, F_overall, &n, &info);
    std::cout << "\033[32mfinish solving matrix;\n\033[0m";
    ;
    // now F_overall is the solution X
    if (info > 0)
    {
        throw Error_throw_ignore("The solution could not be computed!!\n");
    }*/

    /*
    cout << "\nX_overall with BLAS;\n";
    for (size_t i = 0; i < Matrix_D; ++i)
    {
        cout << F_overall[i] << endl;
    }
    */
#pragma omp critical
    {
        std::cout << "\033[32mstart solving matrix;\n\033[0m";
        DFN::Using_UMFPACK U{K_overall, Matrix_D, F_overall};
        std::cout << "\033[32mfinish solving matrix;\n\033[0m";
    }

    this->FEM_results(DFN_mesh, F_overall, dom);

    this->In_and_out_flux(DFN_mesh, dom);
    //cout << "in : " << this->Q_in << endl;
    //cout << "out: " << this->Q_out << endl;

    //delete[] ipiv;
    //ipiv = NULL;
    delete[] K_overall;
    K_overall = NULL;
};

inline void FEM_DFN_coupling::Assemble_overall_matrix(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, double *F_overall, DFN::Domain dom)
{
    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;

    size_t Matrix_D = (NUM_NODES_velocity * 2 + NUM_NODES_p);

    Eigen::MatrixXf KK_big = Eigen::MatrixXf::Zero(Matrix_D, Matrix_D);
    Eigen::VectorXf FF_big = Eigen::VectorXf::Zero(Matrix_D);

    //cout << "NUM_NODES_velocity: " << NUM_NODES_velocity << endl;
    //cout << "NUM_NODES_p: " << NUM_NODES_p << endl;
    //cout << "Matrix_D: " << Matrix_D << endl;
    for (size_t i = 0; i < DFN_mesh.JM_Each_Frac.size(); ++i)
    {
        size_t Frac_Tag = DFN_mesh.Frac_Tag[i];

        double Kper = dom.Fractures[Frac_Tag].Conductivity;
        //cout << "kPer: " << Kper << endl;
        /*
        cout << "Frac " << Frac_Tag << ":\n";
        for(size_t j = 0; j < verts1.size(); ++j)
            cout << verts1[j].transpose() << endl;
        */
        for (size_t j = 0; j < DFN_mesh.JM_Each_Frac[i].size(); ++j)
        {

            Eigen::RowVectorXd w, xi, eta;
            w = Eigen::RowVectorXd::Zero(6);
            w << 0.1099517437,
                0.1099517437,
                0.1099517437,
                0.2233815897,
                0.2233815897,
                0.2233815897;
            xi = Eigen::RowVectorXd::Zero(6);
            xi << 0.0915762135,
                0.0915762135,
                0.8168475730,
                0.4459484909,
                0.4459484909,
                0.1081030182;
            eta = Eigen::RowVectorXd::Zero(6);
            eta << 0.8168475730,
                0.0915762135,
                0.0915762135,
                0.1081030182,
                0.4459484909,
                0.4459484909;

            MatrixXd D1e;
            D1e = MatrixXd::Zero(6, 6);

            Eigen::MatrixXd JXYe_x, JXYe_y, JXYe_x_c, JXYe_y_c;
            JXYe_x = Eigen::MatrixXd::Zero(6, 1);
            JXYe_y = Eigen::MatrixXd::Zero(6, 1);

            for (size_t kk = 0; kk < (size_t)JXYe_x.rows(); ++kk)
            {
                size_t PNT_ID = DFN_mesh.JM_Each_Frac[i][j][kk];
                Vector3d AYUI = DFN_mesh.JXY_2D[i][PNT_ID];

                JXYe_x(kk, 0) = round(AYUI[0], 4);
                JXYe_y(kk, 0) = round(AYUI[1], 4);
            }

            JXYe_x_c = JXYe_x;
            JXYe_y_c = JXYe_y;

            for (size_t ik = 0; ik < (size_t)D1e.rows(); ++ik)
            {
                for (size_t jk = 0; jk < (size_t)D1e.cols(); ++jk)
                {
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

                JXYe_x(kk, 0) = JXYe_x_c(ff, 0);
                JXYe_y(kk, 0) = JXYe_y_c(ff, 0);
            }
            //cout << JXYe_x.transpose() << endl;
            //cout << JXYe_y.transpose() << endl;

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
            //cout << JXYe_x.transpose() << endl;
            //cout << JXYe_y.transpose() << endl;

            for (size_t ik = 0; ik < 6; ++ik)
            {
                for (size_t jk = 0; jk < 6; ++jk)
                {
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
                    size_t Node_m = DFN_mesh.JM_Each_Frac[i][j](jq);
                    size_t Node_n = DFN_mesh.JM_Each_Frac[i][j](kq);

                    KK_big(Node_m, Node_n) += D1e(jq, kq);
                }
            }

            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    size_t Node_m = DFN_mesh.JM_Each_Frac[i][j](jq);
                    size_t Node_n = DFN_mesh.JM_Each_Frac[i][j](kq);

                    KK_big(Node_m + NUM_NODES_velocity, Node_n + NUM_NODES_velocity) += D1e(jq, kq);
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

                    size_t Node_m = DFN_mesh.JM_Each_Frac[i][j](jq);
                    size_t Node_n = DFN_mesh.JM_Each_Frac[i][j](yt);

                    KK_big(Node_m, Node_n + 2 * NUM_NODES_velocity) += C1e(jq, kq);
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

                    size_t Node_m = DFN_mesh.JM_Each_Frac[i][j](jq);
                    size_t Node_n = DFN_mesh.JM_Each_Frac[i][j](yt);

                    /*
                    if (VDF4 == true && Node_m + 1 == 1 && Node_n + 1 == 1)
                        cout << "Node_m: " << Node_m + 1 << ", Node_n: " << Node_n + 1 << endl;*/

                    KK_big(Node_m + NUM_NODES_velocity, Node_n + 2 * NUM_NODES_velocity) += C2e(jq, kq);
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

                    size_t Node_m = DFN_mesh.JM_Each_Frac[i][j](yt);
                    size_t Node_n = DFN_mesh.JM_Each_Frac[i][j](kq);

                    KK_big(Node_m + 2 * NUM_NODES_velocity, Node_n) += B1e(jq, kq);
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

                    size_t Node_m = DFN_mesh.JM_Each_Frac[i][j](yt);
                    size_t Node_n = DFN_mesh.JM_Each_Frac[i][j](kq);

                    KK_big(Node_m + 2 * NUM_NODES_velocity, Node_n + NUM_NODES_velocity) += B2e(jq, kq);
                }
            }
        }
    }

    // boundary condition
    // i = frac
    // j = element
    // 2nd bc
    this->BC_top.resize(DFN_mesh.Frac_Tag.size());
    this->BC_bot.resize(DFN_mesh.Frac_Tag.size());

    for (size_t i = 0; i < DFN_mesh.JM_Each_Frac.size(); ++i)
    {

        for (size_t j = 0; j < DFN_mesh.JM_Each_Frac[i].size(); ++j)
        {
            for (size_t ik = 0; ik < 6; ik += 2)
            {
                size_t edge_0 = ik;
                size_t edge_1 = ik + 1;
                size_t edge_2 = (ik + 2) % 6;

                size_t pnt_ID_0 = DFN_mesh.JM_Each_Frac[i][j](edge_0);
                size_t pnt_ID_1 = DFN_mesh.JM_Each_Frac[i][j](edge_1);
                size_t pnt_ID_2 = DFN_mesh.JM_Each_Frac[i][j](edge_2);

                //-----------find if there are top or bottom edge-------
                if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_top == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_1].If_model_top == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_2].If_model_top == true)
                {
                    BC_top[i].push_back(Vector3d{(double)pnt_ID_0, (double)pnt_ID_1, (double)pnt_ID_2}); //
                    //cout << "1st top: " << pnt_ID_0 + 1 << ", " << pnt_ID_1 + 1 << ", " << pnt_ID_2 + 1 << endl;
                }
                else if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_bottom == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_1].If_model_bottom == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_2].If_model_bottom == true)
                {
                    BC_bot[i].push_back(Vector3d{(double)pnt_ID_0, (double)pnt_ID_1, (double)pnt_ID_2}); //
                }

                // find lateral boundary
                if (((DFN_mesh.Pnt_attri[pnt_ID_0].If_model_top == false || DFN_mesh.Pnt_attri[pnt_ID_1].If_model_top == false || DFN_mesh.Pnt_attri[pnt_ID_2].If_model_top == false) &&
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_bottom == false || DFN_mesh.Pnt_attri[pnt_ID_1].If_model_bottom == false || DFN_mesh.Pnt_attri[pnt_ID_2].If_model_bottom == false)) &&
                    ((DFN_mesh.Pnt_attri[pnt_ID_0].If_model_front == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_front == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_front == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_back == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_back == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_back == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_left == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_left == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_left == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_right == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_right == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_right == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_frac_bound == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_frac_bound == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_frac_bound == true)))
                {

                    //cout << pnt_ID_0 + 1 << ", " << pnt_ID_2 + 1 << endl;
                    double q0 = 0, q1 = 0, q2 = 0; // velocity boundary condition

                    double L = (DFN_mesh.JXY_3D[pnt_ID_0] - DFN_mesh.JXY_3D[pnt_ID_2]).norm();
                    Vector2d Ge;
                    Ge(0) = (L * (q0 + 2 * q1)) / 3;
                    Ge(1) = (L * (2 * q1 + q2)) / 3;

                    size_t Node_s1 = pnt_ID_0, Node_s2 = pnt_ID_2;

                    FF_big[Node_s1 + NUM_NODES_velocity * 2] += Ge(0);
                    FF_big[Node_s2 + NUM_NODES_velocity * 2] += Ge(1);
                    //-------------
                    //cout << "JB2.push_back(Vector4S{" << pnt_ID_0 << "," << pnt_ID_1 << "," << pnt_ID_2 << ", 0});\n";
                }
            }
        }
    }

    // 1st BC
    // i th fracture
    for (size_t i = 0; i < DFN_mesh.NUM_of_linear_NODES; ++i)
    {
        if ((DFN_mesh.Pnt_attri[i].If_model_top == true || DFN_mesh.Pnt_attri[i].If_model_bottom == true))
        {
            double BC_1 = 0;
            if (DFN_mesh.Pnt_attri[i].If_model_top == true)
            {
                BC_1 = 100;
            }

            if (DFN_mesh.Pnt_attri[i].If_model_bottom == true)
            {
                BC_1 = 20;
            }

            size_t node_s = i;
            //cout << "node_s: " << node_s + 1 << ", BC: " << BC_1 << endl;
            //cout << "JB1.push_back(Vector2d(" << node_s << ", " << BC_1 << "));\n";
            node_s += 2 * NUM_NODES_velocity;

            for (size_t go = 0; go < Matrix_D; ++go)
            {
                FF_big[go] -= KK_big(go, node_s) * BC_1;
            }

            for (size_t go = 0; go < Matrix_D; ++go)
            {
                KK_big(go, node_s) = 0; // column
            }

            for (size_t go = 0; go < Matrix_D; ++go)
            {
                KK_big(node_s, go) = 0; // row
            }

            KK_big(node_s, node_s) = 1;
            FF_big[node_s] = BC_1;
        }
    }

    //FF_big = KK_big.colPivHouseholderQr().solve(FF_big);
    //cout << FF_big << endl;
    //exit(0);
    //-----------1st
    for (size_t i = 0; i < Matrix_D; ++i)
    {
        F_overall[i] = FF_big[i];
        for (size_t j = 0; j < Matrix_D; ++j)
        {
            size_t ID_dx = i * Matrix_D + j;
            K_overall[ID_dx] = KK_big(i, j);
        }
    }
}

inline void FEM_DFN_coupling::p_PHI_over_p_x_and_y(double xi,
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

inline void FEM_DFN_coupling::PHI_shape_function(double xi,
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

inline void FEM_DFN_coupling::p_PHI_over_p_xi_and_eta(double xi,
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

inline void FEM_DFN_coupling::p_PSI_over_p_x_and_y(double xi,
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

inline void FEM_DFN_coupling::p_PSI_over_p_xi_and_eta(double xi,
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

inline void FEM_DFN_coupling::PSI_shape_function(double xi,
                                                 double eta,
                                                 double &N_1,
                                                 double &N_2,
                                                 double &N_3)
{
    N_1 = 1 - xi - eta;
    N_2 = xi;
    N_3 = eta;
};

inline void FEM_DFN_coupling::FEM_results(DFN::Mesh_DFN_overall DFN_mesh, double *F_overall, DFN::Domain dom)
{
    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    //
    h_head = Eigen::VectorXd::Zero(NUM_NODES_p);

    Velocity_3D.resize(NUM_NODES_velocity);

    for (size_t i = 0; i < NUM_NODES_p; ++i)
    {
        h_head[i] = F_overall[i + 2 * NUM_NODES_velocity];
    }

    Eigen::VectorXf Index_ve;
    Index_ve = Eigen::VectorXf::Zero(NUM_NODES_velocity);

    for (size_t i = 0; i < DFN_mesh.JM_Each_Frac.size(); ++i)
    {
        size_t FracTAG = DFN_mesh.Frac_Tag[i];
        for (size_t j = 0; j < DFN_mesh.JM_Each_Frac[i].size(); ++j)
        {
            for (size_t k = 0; k < 6; ++k)
            {
                size_t Node_ID = DFN_mesh.JM_Each_Frac[i][j][k];

                if (Index_ve[Node_ID] == 0)
                {
                    Index_ve[Node_ID] = 1;

                    std::vector<Vector3d> v_2d(1), v_3d(1);

                    v_2d[0](0) = F_overall[Node_ID];
                    v_2d[0](1) = F_overall[Node_ID + NUM_NODES_velocity];
                    v_2d[0](2) = 0;

                    DFN::Polygon_convex_3D poly{dom.Fractures[FracTAG].Verts_trim};
                    DFN::Rotate_to_horizontal R1{poly};
                    R1.Rotate_back_without_z(v_2d, v_3d);
                    Velocity_3D[Node_ID] = v_3d[0];
                }
            }
        }
    }
};

inline void FEM_DFN_coupling::In_and_out_flux(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom)
{
    this->Q_in = 0;
    this->Q_out = 0;

    for (size_t i = 0; i < this->BC_top.size(); ++i)
    {
        for (size_t j = 0; j < this->BC_top[i].size(); ++j)
        {
            size_t FRACTAG = DFN_mesh.Frac_Tag[i];

            size_t pointID0 = this->BC_top[i][j](0); // true ID

            //size_t pointID1 = this->BC_top[i][j](1);
            size_t pointID2 = this->BC_top[i][j](2);

            double L = (DFN_mesh.JXY_3D[pointID0] - DFN_mesh.JXY_3D[pointID2]).norm();

            Vector3d V0 = this->Velocity_3D[pointID0];

            Vector3d V2 = this->Velocity_3D[pointID2];

            double ave_velocity = (V0.norm() + V2.norm()) / 2.0;

            double K = dom.Fractures[FRACTAG].Conductivity;

            Q_in += L * K * ave_velocity;
        }
    }

    for (size_t i = 0; i < this->BC_bot.size(); ++i)
    {
        for (size_t j = 0; j < this->BC_bot[i].size(); ++j)
        {
            size_t FRACTAG = DFN_mesh.Frac_Tag[i];

            size_t pointID0 = this->BC_bot[i][j](0); // true ID
            
            //size_t pointID1 = this->BC_bot[i][j](1);
            size_t pointID2 = this->BC_bot[i][j](2);
            
            double L = (DFN_mesh.JXY_3D[pointID0] - DFN_mesh.JXY_3D[pointID2]).norm();

            Vector3d V0 = this->Velocity_3D[pointID0];

            Vector3d V2 = this->Velocity_3D[pointID2];

            double ave_velocity = (V0.norm() + V2.norm()) / 2.0;

            double K = dom.Fractures[FRACTAG].Conductivity;

            Q_out += L * K * ave_velocity;
        }
    }

    cout << "Q_out: " << Q_out << endl;
    cout << "Q_in: " << Q_in << endl;

    if (abs(Q_in - Q_out) / (Q_in > Q_out ? Q_in : Q_out) < 0.3)
    {
        double bc_top = 100, bc_bot = 20;
        double L = dom.Model_domain(0) - dom.Model_domain(1);

        double gradient_head = (bc_top - bc_bot) / L;
        //cout << "gradient_head: " << gradient_head << endl;
        //cout << "Q: " << (Q_out + Q_in) * 0.5 << endl;
        this->Permeability = (Q_out + Q_in) * 0.5 / (pow(L, 2.0) * gradient_head);
        cout << "Permeability " << Permeability << endl;
    }
    else
    {
        cout << "Q_out: " << Q_out << endl;
        cout << "Q_in: " << Q_in << endl;
        //throw Error_throw_ignore("Found large difference between inlet and outlet flux!\n");
    }
};

inline void FEM_DFN_coupling::matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN_overall DFN_mesh, double *F_overall)
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
    //---------------------
    double *pData1;
    double *pData2;
    double *pData3;
    double *pData4;

    pData1 = (double *)mxCalloc(DFN_mesh.JXY_3D.size() * 3, sizeof(double));   // xyz
    pData2 = (double *)mxCalloc(DFN_mesh.JM.size() * 6, sizeof(double));       // topo
    pData3 = (double *)mxCalloc(DFN_mesh.NUM_of_linear_NODES, sizeof(double)); // head
    pData4 = (double *)mxCalloc(DFN_mesh.NUM_of_NODES * 3, sizeof(double));    // velocity_3D

    mxArray *pMxArray1;
    mxArray *pMxArray2;
    mxArray *pMxArray3;
    mxArray *pMxArray4;

    pMxArray1 = mxCreateDoubleMatrix(DFN_mesh.JXY_3D.size(), 3, mxREAL);
    pMxArray2 = mxCreateDoubleMatrix(DFN_mesh.JM.size(), 6, mxREAL);
    pMxArray3 = mxCreateDoubleMatrix(DFN_mesh.NUM_of_linear_NODES, 1, mxREAL);
    pMxArray4 = mxCreateDoubleMatrix(DFN_mesh.NUM_of_NODES, 3, mxREAL);

    if (!pMxArray1 || !pMxArray2 || !pMxArray3 || !pMxArray4)
    {
        throw Error_throw_pause("cannot create pMxArray in class FEM_DFN_coupling\n");
    }

    if (!pData1 || !pData2 || !pData3 || !pData4)
    {
        throw Error_throw_pause("cannot create pData in class FEM_DFN_coupling\n");
    }

    for (size_t j = 0; j < DFN_mesh.JXY_3D.size() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / DFN_mesh.JXY_3D.size()); // column
        l = j % DFN_mesh.JXY_3D.size();       // row

        pData1[j] = DFN_mesh.JXY_3D[l](k);
    }

    for (size_t j = 0; j < DFN_mesh.JM.size() * 6; ++j)
    {
        size_t k, l;
        k = ceil(j / DFN_mesh.JM.size()); // column
        l = j % DFN_mesh.JM.size();       // row

        pData2[j] = DFN_mesh.JM[l](k) + 1;
    }

    for (size_t j = 0; j < DFN_mesh.NUM_of_linear_NODES; ++j)
    {
        pData3[j] = this->h_head[j];
    }

    for (size_t j = 0; j < DFN_mesh.NUM_of_NODES * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / DFN_mesh.NUM_of_NODES); // column
        l = j % DFN_mesh.NUM_of_NODES;       // row

        pData4[j] = this->Velocity_3D[l](k);
    }

    mxSetData(pMxArray1, pData1);
    mxSetData(pMxArray2, pData2);
    mxSetData(pMxArray3, pData3);
    mxSetData(pMxArray4, pData4);

    const char *Frac_JXY_3D = "Frac_JXY_3D";
    const char *Topo_3D = "Topo_3D";
    const char *Head = "Head";
    const char *V_3D = "V_3D";

    matPutVariable(pMatFile, Frac_JXY_3D, pMxArray1);
    matPutVariable(pMatFile, Topo_3D, pMxArray2);
    matPutVariable(pMatFile, Head, pMxArray3);
    matPutVariable(pMatFile, V_3D, pMxArray4);

    mxFree(pData1);
    mxFree(pData2);
    mxFree(pData3);
    mxFree(pData4);
    matClose(pMatFile);

    // m file
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;";
    size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    oss << "%% matrix dimension is " << NUM_NODES_velocity * 2 + NUM_NODES_p << endl;
    oss << "load('" << FileKey_mat << "');\n";
    oss << "figure(1);\n";
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

    oss << "title('DFN head and velocity vector');\n";
    oss << "P = patch('Vertices', " << Frac_JXY_3D << "(1:" << NUM_NODES_p << ",:), 'Faces', " << Topo_3D << "(:, [1,3,5]), 'FaceVertexCData', " << Head << ", 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1);\n";
    oss << "hold on;\n";
    oss << "colorbar;\n";

    oss << "\nquiver3(" << Frac_JXY_3D << "(:,1), " << Frac_JXY_3D << "(:,2), " << Frac_JXY_3D << "(:,3), " << V_3D << "(:,1), " << V_3D << "(:,2), " << V_3D << "(:,3), 'r', 'linewidth', 1.2);\n\n";

    oss.close();
};

inline FEM_DFN_coupling::~FEM_DFN_coupling()
{
    delete[] F_overall;
    F_overall = NULL;
};

}; // namespace DFN
