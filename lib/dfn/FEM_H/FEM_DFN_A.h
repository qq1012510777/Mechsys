#pragma once

//#include "../Mesh_H/Mesh_DFN.h"
#include "../Geometry_H/Normal_vector_2D.h"
#include "../MATLAB_DATA_API/MATLAB_DATA_API.h"
#include "../Mesh_H/Mesh_DFN_overall.h"
#include "../Using_UMFPACK/Using_UMFPACK.h"
namespace DFN
{

class FEM_DFN_A
{
public:
    vector<double> X_overall;
    double Permeability = 0;

    std::vector<Vector3d> Velocity_3D;
    std::vector<Vector2d> Velocity_2D;
    std::vector<Vector3s> BC_IN; // ele ID, node 1, node 2
    std::vector<Vector3s> BC_OUT;
    double Q_in;
    double Q_out;
    double Q_error;

public:
    FEM_DFN_A(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom, size_t i /*direction*/);
    void Assemble_overall_matrix(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, DFN::Domain dom);
    void Apply_boundary_condition(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, double *F_overall, size_t direction, Vector2d BC_head);
    void FEM_results(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom);
    void Identify_In_and_Out_element(DFN::Mesh_DFN_overall DFN_mesh, size_t dir);
    void In_and_Out_flux(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom, size_t dir);
    ~FEM_DFN_A();

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

    void p_PHI_over_p_x_and_y(double xi,
                              double eta,
                              MatrixXd JXYe_x,
                              MatrixXd JXYe_y,
                              VectorXd &pd_N_over_pd_x,
                              VectorXd &pd_N_over_pd_y,
                              MatrixXd &Jacobi);

    void matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN_overall DFN_mesh);
};

inline FEM_DFN_A::FEM_DFN_A(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom, size_t i /*direction*/)
{

    if (i > 2)
    {
        throw Error_throw_pause("Undefined direction!\n");
    }

    size_t Matrix_D = DFN_mesh.NUM_of_NODES;
    Vector2d BC_head;
    BC_head << 100, 20;

    //cout << "\t\t\tMatrix_D = " << Matrix_D << endl;
    //double *K_overall = new double[Matrix_D * Matrix_D]();
    //double *F_overall = new double[Matrix_D]();
    //double K_overall[Matrix_D * Matrix_D] = {};
    //double F_overall[Matrix_D] = {};
    double *K_overall = (double *)malloc(Matrix_D * Matrix_D * sizeof(double));
    double *F_overall = (double *)malloc(Matrix_D * sizeof(double));
    memset(K_overall, 0, Matrix_D * Matrix_D * sizeof(double));
    memset(F_overall, 0, Matrix_D * sizeof(double));
    //cout << "\t\t\tMatrix initialization finished!\n";

    if (K_overall == NULL || F_overall == NULL)
    {
        //delete[] K_overall;
        //K_overall = NULL;
        //delete[] F_overall;
        //F_overall = NULL;
        free(K_overall);
        free(F_overall);
        //cout << "Matrix_D = " << Matrix_D << endl;
        //cout << "Error! Cannot alloc to matrix 'K_overall' or 'F_overall'!\n";
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'K_overall' or 'F_overall'!\n");
    }

    //cout << "\t\t\tassemble matrix start" << endl;
    this->Assemble_overall_matrix(DFN_mesh, K_overall, dom);
    //cout << "\t\t\tassemble matrix finished" << endl;

    //cout << "\t\t\tapply BC start" << endl;
    this->Apply_boundary_condition(DFN_mesh, K_overall, F_overall, i, BC_head);
    //cout << "\t\t\tapply BC finish" << endl;

    //cout << "\t\t\tumfpack start" << endl;
    DFN::Using_UMFPACK U{K_overall, Matrix_D, F_overall};
    //cout << "\t\t\tumfpack finished" << endl;

    X_overall = vector<double>{F_overall, F_overall + Matrix_D};

    //delete[] K_overall;
    //K_overall = NULL;
    //delete[] F_overall;
    //F_overall = NULL;
    free(K_overall);
    free(F_overall);

    this->Identify_In_and_Out_element(DFN_mesh, i);
    this->FEM_results(DFN_mesh, dom);
    this->In_and_Out_flux(DFN_mesh, dom, i);
};

inline FEM_DFN_A::~FEM_DFN_A(){

};

inline void FEM_DFN_A::Assemble_overall_matrix(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, DFN::Domain dom)
{
    size_t Matrix_D = DFN_mesh.NUM_of_NODES;

    for (size_t i = 0; i < DFN_mesh.JM_Each_Frac.size(); ++i)
    {
        size_t Frac_Tag = DFN_mesh.Frac_Tag[i];

        double Kper = dom.Fractures[Frac_Tag].Conductivity;
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

            MatrixXd Ke = Eigen::MatrixXd::Zero(6, 6);

            for (size_t ik = 0; ik < (size_t)Ke.rows(); ++ik)
            {
                for (size_t jk = 0; jk < (size_t)Ke.cols(); ++jk)
                {
                    VectorXd pd_N_over_pd_x, pd_N_over_pd_y;
                    pd_N_over_pd_x = VectorXd::Zero(6);
                    pd_N_over_pd_y = VectorXd::Zero(6);

                    MatrixXd Jacobi;
                    Jacobi = MatrixXd::Zero(2, 2);

                    this->p_PHI_over_p_x_and_y(xi(ik), eta(jk), JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

                    Ke += Kper * w(ik) * w(jk) * Jacobi.determinant() * (pd_N_over_pd_x * pd_N_over_pd_x.transpose() + pd_N_over_pd_y * pd_N_over_pd_y.transpose());
                }
            }

            //---------assemble----------

            for (size_t jq = 0; jq < 6; ++jq)
            {
                for (size_t kq = 0; kq < 6; ++kq)
                {
                    size_t Node_m = DFN_mesh.JM_Each_Frac[i][j](jq);
                    size_t Node_n = DFN_mesh.JM_Each_Frac[i][j](kq);

                    // KK_big(Node_m, Node_n) += Ke(jq, kq);
                    K_overall[Node_m * Matrix_D + Node_n] += Ke(jq, kq);
                }
            }
        }
    }
};

inline void FEM_DFN_A::Apply_boundary_condition(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, double *F_overall, size_t direction, Vector2d BC_head)
{
    size_t Matrix_D = DFN_mesh.NUM_of_NODES;
    //---------------------BC
    /*
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

                // find lateral boundary
                if (((DFN_mesh.Pnt_attri[pnt_ID_0].If_model_top == false || DFN_mesh.Pnt_attri[pnt_ID_1].If_model_top == false || DFN_mesh.Pnt_attri[pnt_ID_2].If_model_top == false) &&
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_bottom == false || DFN_mesh.Pnt_attri[pnt_ID_1].If_model_bottom == false || DFN_mesh.Pnt_attri[pnt_ID_2].If_model_bottom == false)) &&
                    ((DFN_mesh.Pnt_attri[pnt_ID_0].If_model_front == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_front == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_front == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_back == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_back == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_back == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_left == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_left == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_left == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_right == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_model_right == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_model_right == true) ||
                     (DFN_mesh.Pnt_attri[pnt_ID_0].If_frac_bound == true && DFN_mesh.Pnt_attri[pnt_ID_1].If_frac_bound == true && DFN_mesh.Pnt_attri[pnt_ID_2].If_frac_bound == true)))
                {
                    //nothing to do, cause q * n = 0
                }
            }
        }
    }
    */
    for (size_t i = 0; i < DFN_mesh.NUM_of_NODES; ++i)
    {
        double BC_1 = 0;

        bool found_boun_pnt = false;

        if (direction == 2)
        {
            if ((DFN_mesh.Pnt_attri[i].If_model_top == true || DFN_mesh.Pnt_attri[i].If_model_bottom == true))
            {
                found_boun_pnt = true;
                if (DFN_mesh.Pnt_attri[i].If_model_top == true)
                    BC_1 = BC_head[0];

                if (DFN_mesh.Pnt_attri[i].If_model_bottom == true)
                    BC_1 = BC_head[1];
            }
        }
        else if (direction == 1)
        {
            if ((DFN_mesh.Pnt_attri[i].If_model_back == true || DFN_mesh.Pnt_attri[i].If_model_front == true))
            {
                found_boun_pnt = true;
                if (DFN_mesh.Pnt_attri[i].If_model_back == true)
                    BC_1 = BC_head[0];

                if (DFN_mesh.Pnt_attri[i].If_model_front == true)
                    BC_1 = BC_head[1];
            }
        }
        else if (direction == 0)
        {
            if ((DFN_mesh.Pnt_attri[i].If_model_left == true || DFN_mesh.Pnt_attri[i].If_model_right == true))
            {
                found_boun_pnt = true;
                if (DFN_mesh.Pnt_attri[i].If_model_left == true)
                    BC_1 = BC_head[0];

                if (DFN_mesh.Pnt_attri[i].If_model_right == true)
                    BC_1 = BC_head[1];
            }
        }
        else
            throw Error_throw_pause("In FEM_DFN_A, the checking direction is not given!\n");

        if (found_boun_pnt == true)
        {
            int node_s = i;
            for (size_t go = 0; go < Matrix_D; ++go)
            {
                F_overall[go] -= K_overall[go * Matrix_D + node_s] * BC_1;
            }

            for (size_t go = 0; go < Matrix_D; ++go)
            {
                K_overall[go * Matrix_D + node_s] = 0;
            }

            for (size_t go = 0; go < Matrix_D; ++go)
            {
                K_overall[node_s * Matrix_D + go] = 0;
            }

            K_overall[node_s * Matrix_D + node_s] = 1;
            F_overall[node_s] = BC_1;
        }
    }
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

inline void FEM_DFN_A::matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN_overall DFN_mesh)
{
    string Frac_JXY_3D = "Frac_JXY_3D";

    vector<double> JXY_3D_(DFN_mesh.JXY_3D.size() * 3);

    for (size_t j = 0; j < DFN_mesh.JXY_3D.size() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / DFN_mesh.JXY_3D.size()); // column
        l = j % DFN_mesh.JXY_3D.size();       // row

        JXY_3D_[j] = DFN_mesh.JXY_3D[l](k);
    }
    vector<double> JM_(DFN_mesh.JM.size() * 6);
    vector<double> HEAD_(DFN_mesh.NUM_of_NODES);
    vector<double> VELOCITY_(DFN_mesh.JM.size() * 3);

    for (size_t j = 0; j < DFN_mesh.JM.size() * 6; ++j)
    {
        size_t k, l;
        k = ceil(j / DFN_mesh.JM.size()); // column
        l = j % DFN_mesh.JM.size();       // row

        JM_[j] = DFN_mesh.JM[l](k) + 1;
    }

    for (size_t j = 0; j < DFN_mesh.NUM_of_NODES; ++j)
    {
        HEAD_[j] = this->X_overall[j];
    }

    for (size_t j = 0; j < DFN_mesh.JM.size() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / DFN_mesh.JM.size()); // column
        l = j % DFN_mesh.JM.size();       // row

        VELOCITY_[j] = this->Velocity_3D[l](k);
    }

    string Topo_3D = "Topo_3D";
    string Head = "Head";
    string V_3D = "V_3D";

    DFN::MATLAB_DATA_API M1_{FileKey_mat, "w", DFN_mesh.JXY_3D.size() * 3, DFN_mesh.JXY_3D.size(), 3, JXY_3D_, Frac_JXY_3D};
    DFN::MATLAB_DATA_API M2_{FileKey_mat, "u", DFN_mesh.JM.size() * 6, DFN_mesh.JM.size(), 6, JM_, Topo_3D};
    DFN::MATLAB_DATA_API M3_{FileKey_mat, "u", DFN_mesh.NUM_of_NODES, DFN_mesh.NUM_of_NODES, 1, HEAD_, Head};
    DFN::MATLAB_DATA_API M4_{FileKey_mat, "u", DFN_mesh.JM.size() * 3, DFN_mesh.JM.size(), 3, VELOCITY_, V_3D};

    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;";
    //size_t NUM_NODES_velocity = DFN_mesh.NUM_of_NODES;
    //size_t NUM_NODES_p = DFN_mesh.NUM_of_linear_NODES;
    //oss << "%% matrix dimension is " << NUM_NODES_velocity * 2 + NUM_NODES_p << endl;
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

    oss << "title('DFN head');\n";
    oss << "P = patch('Vertices', " << Frac_JXY_3D << ", 'Faces', " << Topo_3D << ", 'FaceVertexCData', " << Head << ", 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1);\n";
    oss << "hold on;\n";
    oss << "colorbar;\n";

    oss << "Center = zeros(" << this->Velocity_3D.size() << ", 3);\n";
    oss << "for i = 1:" << this->Velocity_3D.size() << "\n";
    oss << "\tPnt1 = " << Topo_3D << "(i, 1);\n";
    oss << "\tPnt2 = " << Topo_3D << "(i, 3);\n";
    oss << "\tPnt3 = " << Topo_3D << "(i, 5);\n";
    oss << "\tCenter(i, 1) = (" << Frac_JXY_3D << "(Pnt1, 1) + " << Frac_JXY_3D << "(Pnt2, 1) + " << Frac_JXY_3D << "(Pnt3, 1)) / 3;\n";
    oss << "\tCenter(i, 2) = (" << Frac_JXY_3D << "(Pnt1, 2) + " << Frac_JXY_3D << "(Pnt2, 2) + " << Frac_JXY_3D << "(Pnt3, 2)) / 3;\n";
    oss << "\tCenter(i, 3) = (" << Frac_JXY_3D << "(Pnt1, 3) + " << Frac_JXY_3D << "(Pnt2, 3) + " << Frac_JXY_3D << "(Pnt3, 3)) / 3;\n";
    oss << "end;\n";

    oss << "hold on;\n";
    oss << "quiver3(Center(:,1), Center(:,2), Center(:,3), " << V_3D << "(:,1), " << V_3D << "(:,2), " << V_3D << "(:,3), 'r', 'linewidth', 1.2);\n\n";

    oss.close();
};

inline void FEM_DFN_A::FEM_results(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom)
{
    Velocity_3D.resize(DFN_mesh.JM.size());
    Velocity_2D.resize(DFN_mesh.JM.size());

    /*
    for (size_t i = DFN_mesh.NUM_of_NODES; i < DFN_mesh.NUM_of_NODES + DFN_mesh.JM.size(); ++i)
    {
        //size_t trueID = DFN_mesh.Matrix_PNT_ID[i][2];

        size_t FracTAG = DFN_mesh.JM_Frac_NO[i - DFN_mesh.NUM_of_NODES];
        FracTAG = DFN_mesh.Frac_Tag[FracTAG];

        std::vector<Vector3d> v_2d(1), v_3d(1);

        v_2d[0](0) = this->X_overall[i];
        v_2d[0](1) = this->X_overall[i + DFN_mesh.JM.size()];
        v_2d[0](2) = 0;

        Velocity_2D[i - DFN_mesh.NUM_of_NODES] << v_2d[0](0), v_2d[0](1);

        DFN::Polygon_convex_3D poly{dom.Fractures[FracTAG].Verts_trim};
        DFN::Rotate_to_horizontal R1{poly};
        R1.Rotate_back_without_z(v_2d, v_3d);

        Velocity_3D[i - DFN_mesh.NUM_of_NODES] = v_3d[0];
    }
    */

    for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
    {
        size_t Frac_No = DFN_mesh.JM_Frac_NO[i];

        Eigen::MatrixXd JXYe_x, JXYe_y;
        JXYe_x = Eigen::MatrixXd::Zero(6, 1);
        JXYe_y = Eigen::MatrixXd::Zero(6, 1);

        for (size_t kk = 0; kk < (size_t)JXYe_x.rows(); ++kk)
        {
            JXYe_x(kk, 0) = DFN_mesh.JXY_2D[Frac_No][DFN_mesh.JM[i][kk]](0);
            JXYe_y(kk, 0) = DFN_mesh.JXY_2D[Frac_No][DFN_mesh.JM[i][kk]](1);
        }

        VectorXd pd_N_over_pd_x, pd_N_over_pd_y;
        pd_N_over_pd_x = VectorXd::Zero(6);
        pd_N_over_pd_y = VectorXd::Zero(6);

        MatrixXd Jacobi;
        Jacobi = MatrixXd::Zero(2, 2);

        this->p_PHI_over_p_x_and_y(1. / 3., 1. / 3., JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

        Eigen::VectorXd he;
        he = VectorXd::Zero(6);
        for (size_t j = 0; j < 6; ++j)
            he[j] = this->X_overall[DFN_mesh.JM[i][j]];

        Velocity_2D[i][0] = -1 * dom.Fractures[Frac_No].Conductivity * pd_N_over_pd_x.transpose() * he;
        Velocity_2D[i][1] = -1 * dom.Fractures[Frac_No].Conductivity * pd_N_over_pd_y.transpose() * he;

        //----------------------

        std::vector<Vector3d> v_2d(1), v_3d(1);

        v_2d[0](0) = Velocity_2D[i][0];
        v_2d[0](1) = Velocity_2D[i][1];
        v_2d[0](2) = 0;

        DFN::Polygon_convex_3D poly{dom.Fractures[Frac_No].Verts_trim};
        DFN::Rotate_to_horizontal R1{poly};
        R1.Rotate_back_without_z(v_2d, v_3d);

        Velocity_3D[i] = v_3d[0];
    }
};

inline void FEM_DFN_A::Identify_In_and_Out_element(DFN::Mesh_DFN_overall DFN_mesh, size_t dir)
{
    if (dir == 2)
    {
        for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
        {
            for (size_t ik = 0; ik < 6; ik += 2)
            {
                size_t edge_0 = ik;
                size_t edge_1 = ik + 1;
                size_t edge_2 = (ik + 2) % 6;

                size_t pnt_ID_0 = DFN_mesh.JM[i](edge_0);
                size_t pnt_ID_1 = DFN_mesh.JM[i](edge_1);
                size_t pnt_ID_2 = DFN_mesh.JM[i](edge_2);

                //-----------find if there are top or bottom edge-------
                if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_top == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_1].If_model_top == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_2].If_model_top == true)
                {
                    BC_IN.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
                    //cout << "1st top: " << pnt_ID_0 + 1 << ", " << pnt_ID_1 + 1 << ", " << pnt_ID_2 + 1 << endl;
                }
                else if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_bottom == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_1].If_model_bottom == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_2].If_model_bottom == true)
                {
                    BC_OUT.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
                }
            }
        };
    }
    else if (dir == 1)
    {
        for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
        {
            for (size_t ik = 0; ik < 6; ik += 2)
            {
                size_t edge_0 = ik;
                size_t edge_1 = ik + 1;
                size_t edge_2 = (ik + 2) % 6;

                size_t pnt_ID_0 = DFN_mesh.JM[i](edge_0);
                size_t pnt_ID_1 = DFN_mesh.JM[i](edge_1);
                size_t pnt_ID_2 = DFN_mesh.JM[i](edge_2);

                //-----------find if there are top or bottom edge-------
                if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_back == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_1].If_model_back == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_2].If_model_back == true)
                {
                    BC_IN.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
                    //cout << "1st top: " << pnt_ID_0 + 1 << ", " << pnt_ID_1 + 1 << ", " << pnt_ID_2 + 1 << endl;
                }
                else if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_front == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_1].If_model_front == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_2].If_model_front == true)
                {
                    BC_OUT.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
                }
            }
        };
    }
    else if (dir == 0)
    {
        for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
        {
            for (size_t ik = 0; ik < 6; ik += 2)
            {
                size_t edge_0 = ik;
                size_t edge_1 = ik + 1;
                size_t edge_2 = (ik + 2) % 6;

                size_t pnt_ID_0 = DFN_mesh.JM[i](edge_0);
                size_t pnt_ID_1 = DFN_mesh.JM[i](edge_1);
                size_t pnt_ID_2 = DFN_mesh.JM[i](edge_2);

                //-----------find if there are top or bottom edge-------
                if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_left == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_1].If_model_left == true &&
                    DFN_mesh.Pnt_attri[pnt_ID_2].If_model_left == true)
                {
                    BC_IN.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
                    //cout << "1st top: " << pnt_ID_0 + 1 << ", " << pnt_ID_1 + 1 << ", " << pnt_ID_2 + 1 << endl;
                }
                else if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_right == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_1].If_model_right == true &&
                         DFN_mesh.Pnt_attri[pnt_ID_2].If_model_right == true)
                {
                    BC_OUT.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
                }
            }
        };
    }
    else
        throw Error_throw_pause("Undefined checking direction!\n");
};

inline void FEM_DFN_A::In_and_Out_flux(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom, size_t dir)
{
    Q_in = 0;
    Q_out = 0;
    double L_in = 0, L_out = 0;

    for (size_t i = 0; i < BC_IN.size(); ++i)
    {
        size_t ele_ID = BC_IN[i][0];
        size_t FRACTAG = DFN_mesh.JM_Frac_NO[ele_ID];

        std::vector<Vector2d> lines(2);
        lines[0] << DFN_mesh.JXY_2D[FRACTAG][BC_IN[i][1]][0], DFN_mesh.JXY_2D[FRACTAG][BC_IN[i][1]][1];
        lines[1] << DFN_mesh.JXY_2D[FRACTAG][BC_IN[i][2]][0], DFN_mesh.JXY_2D[FRACTAG][BC_IN[i][2]][1];

        DFN::Normal_vector_2D nor;
        nor.Normal_vector_of_a_line(lines);

        double ave_velocity = abs(this->Velocity_2D[ele_ID].dot(nor.Normal));

        double L = (DFN_mesh.JXY_3D[BC_IN[i][1]] - DFN_mesh.JXY_3D[BC_IN[i][2]]).norm();

        double K = dom.Fractures[FRACTAG].Conductivity;

        Q_in += L * K * ave_velocity;

        L_in += L;
    }
    //cout << "---------------------\n";
    for (size_t i = 0; i < BC_OUT.size(); ++i)
    {
        size_t ele_ID = BC_OUT[i][0];
        size_t FRACTAG = DFN_mesh.JM_Frac_NO[ele_ID];

        std::vector<Vector2d> lines(2);
        lines[0] << DFN_mesh.JXY_2D[FRACTAG][BC_OUT[i][1]][0], DFN_mesh.JXY_2D[FRACTAG][BC_OUT[i][1]][1];
        lines[1] << DFN_mesh.JXY_2D[FRACTAG][BC_OUT[i][2]][0], DFN_mesh.JXY_2D[FRACTAG][BC_OUT[i][2]][1];
        DFN::Normal_vector_2D nor;
        nor.Normal_vector_of_a_line(lines);

        double ave_velocity = abs(this->Velocity_2D[ele_ID].dot(nor.Normal));

        double L = (DFN_mesh.JXY_3D[BC_OUT[i][1]] - DFN_mesh.JXY_3D[BC_OUT[i][2]]).norm();

        double K = dom.Fractures[FRACTAG].Conductivity;

        Q_out += L * K * ave_velocity;
        L_out += L;
    }

    //cout << "Difference: " << (abs(Q_in - Q_out) / (Q_in > Q_out ? Q_in : Q_out)) * 100 << "%" << endl;
    //cout << "Q_in: " << Q_in << ", Q_out: " << Q_out << endl;
    //cout << "L_in: " << L_in << ", L_out: " << L_out << endl;

    /*
    if ((abs(Q_in - Q_out) / (Q_in > Q_out ? Q_in : Q_out)) > 0.1)
    {
        string AS = "Warning! The difference between Q_in and Q_out is too large!\n";
        AS = AS + "Q_in: " + to_string(Q_in) + "\n";
        AS = AS + "Q_out: " + to_string(Q_out) + "\n";
        AS = AS + "Difference: " + to_string(abs(Q_in - Q_out) / (Q_in > Q_out ? Q_in : Q_out) * 100) + " %\n";
        throw Error_throw_ignore(AS);
    }*/

    this->Q_error = (abs(Q_in - Q_out) / (Q_in > Q_out ? Q_in : Q_out)) * 100.;

    double l = 0, m = 0;
    if (dir == 0)
    {
        l = 4;
        m = 5;
    }
    else if (dir == 1)
    {
        l = 2;
        m = 3;
    }
    else if (dir == 2)
    {
        l = 1;
        m = 0;
    }

    this->Permeability = 0.5 * (Q_out + Q_in) / abs(dom.Model_domain[l] - dom.Model_domain[m]);
}
}; // namespace DFN
