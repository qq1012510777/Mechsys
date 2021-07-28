#pragma once

//#include "../Mesh_H/Mesh_DFN.h"
#include "../Geometry_H/Normal_vector_2D.h"
#include "../Mesh_H/Mesh_DFN_overall.h"
#include "../Using_UMFPACK/Using_UMFPACK.h"
namespace DFN
{

class FEM_DFN_A
{
public:
    double *F_overall;
    double Permeability = 0;
    std::vector<Vector3d> Velocity_3D;
    std::vector<Vector2d> Velocity_2D;
    std::vector<Vector3s> BC_TOP; // ele ID, node 1, node 2
    std::vector<Vector3s> BC_BOT;
    double Q_in;
    double Q_out;

public:
    FEM_DFN_A(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom);
    void Assemble_overall_matrix(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, double *F_overall, DFN::Domain dom);
    void FEM_results(DFN::Mesh_DFN_overall DFN_mesh, double *F_overall, DFN::Domain dom);
    void In_and_Out_flux(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom);
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

    void matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN_overall DFN_mesh, double *F_overall);
};

inline FEM_DFN_A::FEM_DFN_A(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom)
{
    size_t Matrix_D = DFN_mesh.NUM_of_NODES;
    //size_t NUM_ELES = DFN_mesh.JM.size();
    //size_t NUM_TRACE_ELES = DFN_mesh.NUM_trace_ele_sets;

    //Matrix_D += (2 * NUM_ELES /*+ NUM_TRACE_ELES*/);

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

#pragma omp critical
    {
        std::cout << "\033[32mstart solving matrix;\n\033[0m";
        DFN::Using_UMFPACK U{K_overall, Matrix_D, F_overall};
        std::cout << "\033[32mfinish solving matrix;\n\033[0m";
    }

    /*  
    cout << "\nF_overall;\n";
    for (size_t i = 0; i < Matrix_D; ++i)
    {
        cout << F_overall[i] << endl;
    }
    */
    this->FEM_results(DFN_mesh, F_overall, dom);
    this->In_and_Out_flux(DFN_mesh, dom);
    delete[] K_overall;
    K_overall = NULL;
};

inline FEM_DFN_A::~FEM_DFN_A()
{
    delete[] F_overall;
    F_overall = NULL;
};

inline void FEM_DFN_A::Assemble_overall_matrix(DFN::Mesh_DFN_overall DFN_mesh, double *K_overall, double *F_overall, DFN::Domain dom)
{
    size_t Matrix_D = DFN_mesh.NUM_of_NODES;
    //size_t NUM_ELES = DFN_mesh.JM.size();
    //size_t NUM_TRACE_ELES = DFN_mesh.NUM_trace_ele_sets;

    //Matrix_D += (2 * NUM_ELES /*+ NUM_TRACE_ELES*/);

    Eigen::MatrixXf KK_big = Eigen::MatrixXf::Zero(Matrix_D, Matrix_D);
    Eigen::VectorXf FF_big = Eigen::VectorXf::Zero(Matrix_D);

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

                    KK_big(Node_m, Node_n) += Ke(jq, kq);
                }
            }
        }
    }

    //---------------------BC
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

    for (size_t i = 0; i < DFN_mesh.NUM_of_NODES; ++i)
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

            int node_s = i;
            //cout << "node_s: " << node_s + 1 << ", BC: " << BC_1 << endl;
            //cout << "JB1.push_back(Vector2d(" << node_s << ", " << BC_1 << "));\n";

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

    for (size_t i = 0; i < Matrix_D; ++i)
    {
        F_overall[i] = FF_big[i];
        for (size_t j = 0; j < Matrix_D; ++j)
        {
            size_t ID_dx = i * Matrix_D + j;
            K_overall[ID_dx] = KK_big(i, j);
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

inline void FEM_DFN_A::matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom, DFN::Mesh_DFN_overall DFN_mesh, double *F_overall)
{

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

    pData1 = (double *)mxCalloc(DFN_mesh.JXY_3D.size() * 3, sizeof(double)); // xyz
    pData2 = (double *)mxCalloc(DFN_mesh.JM.size() * 6, sizeof(double));     // topo
    pData3 = (double *)mxCalloc(DFN_mesh.NUM_of_NODES, sizeof(double));      // head
    pData4 = (double *)mxCalloc(DFN_mesh.JM.size() * 3, sizeof(double));     // head

    mxArray *pMxArray1;
    mxArray *pMxArray2;
    mxArray *pMxArray3;
    mxArray *pMxArray4;

    pMxArray1 = mxCreateDoubleMatrix(DFN_mesh.JXY_3D.size(), 3, mxREAL);
    pMxArray2 = mxCreateDoubleMatrix(DFN_mesh.JM.size(), 6, mxREAL);
    pMxArray3 = mxCreateDoubleMatrix(DFN_mesh.NUM_of_NODES, 1, mxREAL);
    pMxArray4 = mxCreateDoubleMatrix(DFN_mesh.JM.size(), 3, mxREAL);

    if (!pMxArray1 || !pMxArray2 || !pMxArray3 || !pMxArray4)
    {
        throw Error_throw_pause("cannot create pMxArray in class FEM_DFN_A\n");
    }

    if (!pData1 || !pData2 || !pData3 || !pData4)
    {
        throw Error_throw_pause("cannot create pData in class FEM_DFN_A\n");
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

    for (size_t j = 0; j < DFN_mesh.NUM_of_NODES; ++j)
    {
        pData3[j] = this->F_overall[j];
    }

    for (size_t j = 0; j < DFN_mesh.JM.size() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / DFN_mesh.JM.size()); // column
        l = j % DFN_mesh.JM.size();       // row

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

inline void FEM_DFN_A::FEM_results(DFN::Mesh_DFN_overall DFN_mesh, double *F_overall, DFN::Domain dom)
{
    Velocity_3D.resize(DFN_mesh.JM.size());
    Velocity_2D.resize(DFN_mesh.JM.size());

    for (size_t i = DFN_mesh.NUM_of_NODES; i < DFN_mesh.NUM_of_NODES + DFN_mesh.JM.size(); ++i)
    {
        //size_t trueID = DFN_mesh.Matrix_PNT_ID[i][2];

        size_t FracTAG = DFN_mesh.JM_Frac_NO[i - DFN_mesh.NUM_of_NODES];
        FracTAG = DFN_mesh.Frac_Tag[FracTAG];

        std::vector<Vector3d> v_2d(1), v_3d(1);

        v_2d[0](0) = F_overall[i];
        v_2d[0](1) = F_overall[i + DFN_mesh.JM.size()];
        v_2d[0](2) = 0;

        Velocity_2D[i - DFN_mesh.NUM_of_NODES] << v_2d[0](0), v_2d[0](1);

        DFN::Polygon_convex_3D poly{dom.Fractures[FracTAG].Verts_trim};
        DFN::Rotate_to_horizontal R1{poly};
        R1.Rotate_back_without_z(v_2d, v_3d);

        Velocity_3D[i - DFN_mesh.NUM_of_NODES] = v_3d[0];
    }

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
            he[j] = F_overall[DFN_mesh.JM[i][j]];

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

inline void FEM_DFN_A::In_and_Out_flux(DFN::Mesh_DFN_overall DFN_mesh, DFN::Domain dom)
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
                BC_TOP.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
                //cout << "1st top: " << pnt_ID_0 + 1 << ", " << pnt_ID_1 + 1 << ", " << pnt_ID_2 + 1 << endl;
            }
            else if (DFN_mesh.Pnt_attri[pnt_ID_0].If_model_bottom == true &&
                     DFN_mesh.Pnt_attri[pnt_ID_1].If_model_bottom == true &&
                     DFN_mesh.Pnt_attri[pnt_ID_2].If_model_bottom == true)
            {
                BC_BOT.push_back(Vector3s{i, pnt_ID_0, pnt_ID_2}); //
            }
        }
    };
    Q_in = 0;
    Q_out = 0;
    double L_in = 0, L_out = 0;

    for (size_t i = 0; i < BC_TOP.size(); ++i)
    {
        size_t ele_ID = BC_TOP[i][0];
        size_t FRACTAG = DFN_mesh.JM_Frac_NO[ele_ID];

        std::vector<Vector2d> lines(2);
        lines[0] << DFN_mesh.JXY_2D[FRACTAG][BC_TOP[i][1]][0], DFN_mesh.JXY_2D[FRACTAG][BC_TOP[i][1]][1];
        lines[1] << DFN_mesh.JXY_2D[FRACTAG][BC_TOP[i][2]][0], DFN_mesh.JXY_2D[FRACTAG][BC_TOP[i][2]][1];

        DFN::Normal_vector_2D nor;
        nor.Normal_vector_of_a_line(lines);

        double ave_velocity = abs(this->Velocity_2D[ele_ID].dot(nor.Normal));

        double L = (DFN_mesh.JXY_3D[BC_TOP[i][1]] - DFN_mesh.JXY_3D[BC_TOP[i][2]]).norm();

        double K = dom.Fractures[FRACTAG].Conductivity;

        Q_in += L * K * ave_velocity;

        L_in += L;
    }
    //cout << "---------------------\n";
    for (size_t i = 0; i < BC_BOT.size(); ++i)
    {
        size_t ele_ID = BC_BOT[i][0];
        size_t FRACTAG = DFN_mesh.JM_Frac_NO[ele_ID];

        std::vector<Vector2d> lines(2);
        lines[0] << DFN_mesh.JXY_2D[FRACTAG][BC_BOT[i][1]][0], DFN_mesh.JXY_2D[FRACTAG][BC_BOT[i][1]][1];
        lines[1] << DFN_mesh.JXY_2D[FRACTAG][BC_BOT[i][2]][0], DFN_mesh.JXY_2D[FRACTAG][BC_BOT[i][2]][1];
        DFN::Normal_vector_2D nor;
        nor.Normal_vector_of_a_line(lines);

        double ave_velocity = abs(this->Velocity_2D[ele_ID].dot(nor.Normal));

        double L = (DFN_mesh.JXY_3D[BC_BOT[i][1]] - DFN_mesh.JXY_3D[BC_BOT[i][2]]).norm();

        double K = dom.Fractures[FRACTAG].Conductivity;

        Q_out += L * K * ave_velocity;
        L_out += L;
    }

    cout << "Q_in: " << Q_in << endl;
    cout << "Q_out: " << Q_out << endl
         << endl;

    cout << "L_in: " << L_in << endl;
    cout << "L_out: " << L_out << endl;
}

}; // namespace DFN
