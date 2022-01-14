#pragma once
#include "../Error_throw/Error_throw.h"
#include "../Geometry_H/Normal_vector_2D.h"
#include "../MATLAB_DATA_API/MATLAB_DATA_API.h"

#include "../Mesh_H/Mesh_DFN_linear.h"
//#include "../Using_UMFPACK/Using_UMFPACK.h"
//#include "Eigen/CholmodSupport"
//#include <Eigen/UmfPackSupport>

namespace DFN
{
class MHFEM
{
public:
    double inlet_p = 0;
    double outlet_p = 0;

    double Permeability = 0;
    double Q_in = 0;
    double Q_out = 0;
    double Q_error = 0;

    double inlet_length = 0;
    double outlet_length = 0;
    size_t dir = 0;

    MatrixXd velocity_normal_scalar_sep_edges;
    MatrixXd pressure_interior_edge;
    MatrixXd pressure_eles;

    MatrixXf normal_vec_sep_edges;

    size_t round_precision = 10;

    size_t N_proc = 1;

public:
    MHFEM(DFN::Mesh_DFN_linear mesh,
          DFN::Domain dom,
          double inlet_p_,
          double outlet_p_,
          size_t i /*direction*/,
          size_t no_proc = 1);

    void Matlab_plot(string FileKey_mat,
                     string FileKey_m,
                     DFN::Mesh_DFN_linear mesh,
                     DFN::Domain dom,
                     string same_dir = "NO");

private:
    void Assemble_and_solve(DFN::Mesh_DFN_linear mesh,
                            DFN::Domain dom);

    MatrixXf Calculate_velocity_normal_vec_sep_edges(DFN::Domain dom,
                                                     size_t Frac_ID,
                                                     MatrixXf coord);

    double The_inlet_of_all_elements(DFN::Mesh_DFN_linear mesh); // do not use

private:
    MatrixXf stimaB(MatrixXf coord);
    size_t Element_ID(size_t i, size_t j, DFN::Mesh_DFN_linear mesh);
};

inline MHFEM::MHFEM(DFN::Mesh_DFN_linear mesh,
                    DFN::Domain dom,
                    double inlet_p_,
                    double outlet_p_,
                    size_t i /*direction*/,
                    size_t no_proc)
{
    this->round_precision = 10;

    dir = i;
    inlet_p = inlet_p_;
    outlet_p = outlet_p_;
    this->N_proc = no_proc;
    Assemble_and_solve(mesh,
                       dom);

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
};

inline void MHFEM::Assemble_and_solve(DFN::Mesh_DFN_linear mesh,
                                      DFN::Domain dom)
{
    //auto start = std::chrono::steady_clock::now(), end = std::chrono::steady_clock::now();

    size_t NUM_sep_edges = mesh.element_3D.rows() * 3,
           NUM_eles = mesh.element_3D.rows(),
           NUM_glob_interior_edges = mesh.NUM_interior_edges;

    size_t Dim = NUM_sep_edges + NUM_eles + NUM_glob_interior_edges;

    normal_vec_sep_edges.resize(NUM_sep_edges, 3);

    SparseMatrix<double> K(Dim, Dim);
    SparseMatrix<double> b(Dim, 1);

    K.reserve(VectorXi::Constant(Dim, 5));
    b.reserve(mesh.NUM_inlet_edges + mesh.NUM_outlet_edges);

    vector<Triplet<double>> Inlet_edge_sep_NO,
        Outlet_edge_sep_NO;

    vector<double> length_in_out(mesh.NUM_inlet_edges + mesh.NUM_outlet_edges, 0);

    vector<size_t> Neumann_edge_sep_NO;

    Inlet_edge_sep_NO.resize(mesh.NUM_inlet_edges);
    Outlet_edge_sep_NO.resize(mesh.NUM_outlet_edges);
    Neumann_edge_sep_NO.resize(mesh.NUM_neumann_edges);

    //cout << "*********assembling_1*********\n";
    //time_counter_start(start);
    size_t Nproc_1 = this->N_proc;

    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i) // Frac
    {
        size_t Frac_Tag = mesh.Frac_Tag[i];

        double Kper = dom.Fractures[Frac_Tag].Conductivity;

#pragma omp parallel for schedule(dynamic) num_threads(Nproc_1)
        for (size_t j = 0; j < (size_t)mesh.element_2D[i].rows(); ++j) // elements
        {

            size_t tmp_eleNO = Element_ID(i, j, mesh);

            //cout << tmp_eleNO << endl;

            MatrixXf coord = MatrixXf::Zero(3, 2);

            size_t node1 = mesh.element_2D[i](j, 0),
                   node2 = mesh.element_2D[i](j, 1),
                   node3 = mesh.element_2D[i](j, 2);

            coord.row(0) << mesh.coordinate_2D[i].coeffRef(node1 - 1, 0),
                mesh.coordinate_2D[i].coeffRef(node1 - 1, 1);

            coord.row(1) << mesh.coordinate_2D[i].coeffRef(node2 - 1, 0),
                mesh.coordinate_2D[i].coeffRef(node2 - 1, 1);

            coord.row(2) << mesh.coordinate_2D[i].coeffRef(node3 - 1, 0),
                mesh.coordinate_2D[i].coeffRef(node3 - 1, 1);

            MatrixXf A_ = stimaB(coord);

            size_t edge1 = tmp_eleNO * 3 + 1,
                   edge2 = tmp_eleNO * 3 + 2,
                   edge3 = tmp_eleNO * 3 + 3;

            /*
            #pragma omp critical
            if (i == 32)
            {
                cout << "ele: " << endl;
                cout << tmp_eleNO + 1 << endl;
                cout << "sep velovity:" << endl;
                cout << edge1 << endl;
                cout << edge2 << endl;
                cout << edge3 << endl;
            }*/

            MatrixXf outnormal_vec = Calculate_velocity_normal_vec_sep_edges(dom, Frac_Tag, coord);
            normal_vec_sep_edges.block(edge1 - 1, 0, 3, 3) = outnormal_vec;

            //vector<size_t> I{edge2 - 1, edge3 - 1, edge1 - 1};
            vector<size_t> I{edge2 - 1, edge3 - 1, edge1 - 1};

            float SUM_u = 0;
            for (size_t ik = 0; ik < 3; ++ik)
                for (size_t jk = 0; jk < 3; ++jk)
                {
                    double AT = (double)round(((-1.0 / Kper) * A_(ik, jk)), round_precision);
                    //#pragma omp critical
                    {
                        K.insert(I[ik], I[jk]) = AT;
                    }
                    SUM_u += abs(AT);
                }

            //cout << abs(SUM_u) << endl;
            if (abs(SUM_u) < pow(10, -1.0 * round_precision))
            {
                string AS = "bad triangle lead to all zero matrix A_!\n";
                AS += "Frac: " + to_string(i + 1) + "; element: " + to_string(j + 1) + "\n\nThe matrix is:\n";
                vector<RowVector3f> coord_(3);
                coord_[0] = mesh.coordinate_3D.row(node1 - 1);
                coord_[1] = mesh.coordinate_3D.row(node2 - 1);
                coord_[2] = mesh.coordinate_3D.row(node3 - 1);

                DFN::If_skinny_triangle SKINNY{coord_};

                for (size_t io = 0; io < 3; ++io)
                {
                    for (size_t jo = 0; jo < 3; ++jo)
                    {
                        AS += to_string(A_(io, jo));
                        AS += ", ";
                    }
                    AS += "\n";
                }

                AS += "\nAngle is:\n";
                AS += to_string(SKINNY.angle);
                AS += "\n";
                AS += "abs(SUM_u): ";
                AS += to_string(abs(SUM_u));
                AS += "\n";
                cout << AS << endl;
                throw Error_throw_ignore(AS);
            }

            Vector3f B_;
            B_ << round((coord.row(2) - coord.row(1)).norm(), round_precision),
                round((coord.row(0) - coord.row(2)).norm(), round_precision),
                round((coord.row(1) - coord.row(0)).norm(), round_precision);

            /*
            if (i == 16)
            {
            #pragma omp critical
                {
                    cout << "element: " << j + 1 << endl;
                    cout << A_ << endl
                         << endl;
                    cout << B_ << endl;
                }
            }*/

            for (size_t ik = 0; ik < 3; ++ik)
            {
                //B.coeffRef(I[ik], tmp_eleNO) = (double)round(B_(ik, 0), round_precision);

                //#pragma omp critical
                {
                    K.insert(I[ik], tmp_eleNO + NUM_sep_edges) = (double)round(B_(ik, 0), round_precision);
                    K.insert(tmp_eleNO + NUM_sep_edges, I[ik]) = (double)round(B_(ik, 0), round_precision);
                }

                if (mesh.Interior_edgeNO[i](j, ik).first == "interior")
                {
                    size_t global_interiorID = mesh.Interior_edgeNO[i](j, ik).second[0];

                    K.insert(global_interiorID - 1 + NUM_sep_edges + NUM_eles, I[(ik + 2) % 3]) = -B_((ik + 2) % 3, 0);
                    K.insert(I[(ik + 2) % 3], global_interiorID - 1 + NUM_sep_edges + NUM_eles) = -B_((ik + 2) % 3, 0);
                    /*
                    #pragma omp critical
                    if (i == 32)
                    {
                        cout << "interior: " << endl;
                        cout << global_interiorID << endl;
                    }*/
                }
                else if (mesh.Interior_edgeNO[i](j, ik).first == "in")
                {
                    //size_t node_a = mesh.element_2D[i](j, ik);
                    //size_t node_b = mesh.element_2D[i](j, (ik + 1) % 3);

                    size_t Sep_NO_rr = mesh.Interior_edgeNO[i](j, ik).second[0];
                    size_t ID_1 = mesh.Interior_edgeNO[i](j, ik).second[1];

                    //Inlet_edge_sep_NO[ID_1 - 1].first = Sep_NO_rr;
                    //Inlet_edge_sep_NO[ID_1 - 1].second = (double)round(B_((ik + 2) % 3, 0), round_precision);
                    double P_D = this->inlet_p;

                    Inlet_edge_sep_NO[ID_1 - 1] = Triplet<double>(Sep_NO_rr - 1, 0, P_D * (double)round(B_((ik + 2) % 3, 0), round_precision));
                    length_in_out[ID_1 - 1] = (double)round(B_((ik + 2) % 3, 0), round_precision);
                }
                else if (mesh.Interior_edgeNO[i](j, ik).first == "out")
                {
                    //size_t node_a = mesh.element_2D[i](j, ik);
                    //size_t node_b = mesh.element_2D[i](j, (ik + 1) % 3);

                    size_t Sep_NO_rr = mesh.Interior_edgeNO[i](j, ik).second[0];
                    size_t ID_1 = mesh.Interior_edgeNO[i](j, ik).second[1];

                    // Outlet_edge_sep_NO[ID_1 - 1].first = Sep_NO_rr;
                    // Outlet_edge_sep_NO[ID_1 - 1].second = (double)round(B_((ik + 2) % 3, 0), round_precision);

                    double P_D = this->outlet_p;
                    Outlet_edge_sep_NO[ID_1 - 1] = Triplet<double>(Sep_NO_rr - 1, 0, P_D * (double)round(B_((ik + 2) % 3, 0), round_precision));
                    length_in_out[ID_1 - 1 + mesh.NUM_inlet_edges] = (double)round(B_((ik + 2) % 3, 0), round_precision);
                }
                else if (mesh.Interior_edgeNO[i](j, ik).first == "neumann")
                {
                    //size_t node_a = mesh.element_2D[i](j, ik);
                    //size_t node_b = mesh.element_2D[i](j, (ik + 1) % 3);

                    size_t Sep_NO_rr = mesh.Interior_edgeNO[i](j, ik).second[0];
                    size_t ID_1 = mesh.Interior_edgeNO[i](j, ik).second[1];

                    Neumann_edge_sep_NO[ID_1 - 1] = Sep_NO_rr;
                    //Neumann_edge_sep_NO[ID_1 - 1].second = (double)round(B_((ik + 2) % 3, 0), round_precision);
                }
            }
        }
    }
    //time_counter_end(start, end, "assembling_1", "in_minutes");

    //cout << "*********boundary condition*********\n";
    //time_counter_start(start);
    Eigen::setNbThreads(this->N_proc);
    vector<Triplet<double>> tripletlist;
    tripletlist.reserve(Inlet_edge_sep_NO.size() + Outlet_edge_sep_NO.size());
    tripletlist.insert(tripletlist.end(), Inlet_edge_sep_NO.begin(), Inlet_edge_sep_NO.end());
    tripletlist.insert(tripletlist.end(), Outlet_edge_sep_NO.begin(), Outlet_edge_sep_NO.end());

    b.setFromTriplets(tripletlist.begin(), tripletlist.end());

    //--------Neumann
    for (size_t i = 0; i < Neumann_edge_sep_NO.size(); ++i)
    {
        size_t SEP_edgeNO = Neumann_edge_sep_NO[i];
        K.innerVector(SEP_edgeNO - 1) = VectorXd::Zero(Dim, 1).sparseView();

        K.row(SEP_edgeNO - 1) *= 0;

        K.coeffRef(SEP_edgeNO - 1, SEP_edgeNO - 1) = 1;
        //b.coeffRef(SEP_edgeNO - 1, 0) = 0; // impermeable
    }

    K.makeCompressed();
    b.makeCompressed();

    SparseMatrix<double> A(NUM_sep_edges, NUM_sep_edges);
    SparseMatrix<double> B(NUM_eles, NUM_sep_edges);
    SparseMatrix<double> C(NUM_glob_interior_edges, NUM_sep_edges);

    A = K.block(0, 0, NUM_sep_edges, NUM_sep_edges);
    B = K.block(NUM_sep_edges, 0, NUM_eles, NUM_sep_edges);
    C = K.block(NUM_sep_edges + NUM_eles, 0, NUM_glob_interior_edges, NUM_sep_edges);

    K.resize(0, 0);

    A.makeCompressed();
    B.makeCompressed();
    C.makeCompressed();

    SparseMatrix<double> g = b.block(0, 0, NUM_sep_edges, 1);

    SparseMatrix<double> f = b.block(NUM_sep_edges, 0, NUM_eles, 1);

    b.resize(0, 0);

    g.makeCompressed();
    f.makeCompressed();
    //time_counter_end(start, end, "boundary condition addressing", "in_minutes");

    //cout << "*********preparing*********\n";
    //time_counter_start(start);
    SparseMatrix<double> A_inv(A.rows(), A.cols());
    A_inv.reserve(VectorXi::Constant(A.cols(), 3));

#pragma omp parallel for schedule(dynamic) num_threads(Nproc_1)
    for (int i = 0; i < A.rows() / 3; ++i)
    {
        MatrixXd A_block = MatrixXd(A.block(i * 3, i * 3, 3, 3));
        //A_inv.block(i * 3, i * 3, 3, 3) = A_block.inverse();
        A_block = A_block.inverse();
        for (int j = i * 3, jq = 0; j < i * 3 + 3; ++j, ++jq)
            for (int k = i * 3, kq = 0; k < i * 3 + 3; ++k, ++kq)
                A_inv.insert(j, k) = A_block(jq, kq);
    }
    A_inv.makeCompressed();
    A.resize(0, 0);

    SparseMatrix<double> C_tps = C.transpose();
    C_tps.makeCompressed();
    SparseMatrix<double> B_tps = B.transpose();
    B_tps.makeCompressed();
    SparseMatrix<double> Wq = B * A_inv;
    Wq.makeCompressed();
    B.resize(0, 0);

    SparseMatrix<double> U = Wq * B_tps;

    //cout << 1 << endl;
#pragma omp parallel for schedule(dynamic) num_threads(Nproc_1)
    for (int i = 0; i < U.rows(); ++i)
        U.coeffRef(i, i) = 1.0 / U.coeffRef(i, i);
    //cout << 1 << endl;
    U.makeCompressed();

    SparseMatrix<double> Re = C * A_inv;
    Re.makeCompressed();
    C.resize(0, 0);
    SparseMatrix<double> Eq = Re * B_tps;
    Eq.makeCompressed();
    SparseMatrix<double> Sd = Wq * C_tps;
    Sd.makeCompressed();

    SparseMatrix<double> D = Re * C_tps - Eq * U * Sd;
    D.makeCompressed();
    SparseMatrix<double> r = Re * g + Eq * U * (f - Wq * g);
    r.makeCompressed();
    //time_counter_end(start, end, "preparing matrice", "in_minutes");

    //cout << "*********solving*********\n";
    //time_counter_start(start);
    SimplicialLDLT<SparseMatrix<double>> solver;
    pressure_interior_edge = solver.compute(D).solve(r);
    // time_counter_end(start, end, "solving matrice", "in_minutes");

    pressure_eles = U * (Wq * g - Sd * pressure_interior_edge - f);
    velocity_normal_scalar_sep_edges = A_inv * (g - B_tps * pressure_eles - C_tps * pressure_interior_edge);

    //cout << "in:\n";
    for (size_t i = 0; i < Inlet_edge_sep_NO.size(); ++i)
    {
        size_t sep_EDGE_no = Inlet_edge_sep_NO[i].row();
        //double len = Inlet_edge_sep_NO[i].value();
        double len = length_in_out[i];

        double veloc_length = abs(velocity_normal_scalar_sep_edges(sep_EDGE_no, 0) * len);
        this->Q_in += veloc_length;
        inlet_length += len;
        //cout << "len: " << len << ";\tq: " << velocity_normal_scalar_sep_edges(sep_EDGE_no, 0) << "; sep_EDGE_no: " << sep_EDGE_no + 1 << "\n";
    }
    //cout << "\n\nout:\n";
    for (size_t i = 0; i < Outlet_edge_sep_NO.size(); ++i)
    {
        size_t sep_EDGE_no = Outlet_edge_sep_NO[i].row();
        //double len = Outlet_edge_sep_NO[i].value();
        double len = length_in_out[i + mesh.NUM_inlet_edges];

        double veloc_length = abs(velocity_normal_scalar_sep_edges(sep_EDGE_no, 0) * len);
        this->Q_out += veloc_length;
        outlet_length += len;
        //cout << "len: " << len << ";\tq: " << velocity_normal_scalar_sep_edges(sep_EDGE_no, 0) << "; sep_EDGE_no: " << sep_EDGE_no + 1 << "\n";
    }

}; // namespace DFN

inline MatrixXf MHFEM::stimaB(MatrixXf coord)
{
    Matrix<float, Dynamic, Dynamic> N;
    N = MatrixXf::Zero(6, 3);

    VectorXf P1_P2, P2_P1, P3_P2, P2_P3, P3_P1, P1_P3;
    P1_P2 = coord.row(0) - coord.row(1);
    P2_P1 = -P1_P2;
    P3_P2 = coord.row(2) - coord.row(1);
    P2_P3 = -P3_P2;
    P3_P1 = coord.row(2) - coord.row(0);
    P1_P3 = -P3_P1;

    N.col(0) << 0, 0, P2_P1[0], P2_P1[1], P3_P1[0], P3_P1[1];
    N.col(1) << P1_P2[0], P1_P2[1], 0, 0, P3_P2[0], P3_P2[1];
    N.col(2) << P1_P3[0], P1_P3[1], P2_P3[0], P2_P3[1], 0, 0;

    MatrixXf D = MatrixXf::Zero(3, 3);
    D(0, 0) = P3_P2.norm();
    D(1, 1) = P1_P3.norm();
    D(2, 2) = P1_P2.norm();

    MatrixXf M, M1, M2;
    M = MatrixXf::Zero(6, 6);

    M1 = MatrixXf::Zero(3, 3);
    M2 = M1;

    M1.row(0) << 2, 0, 1;
    M1.row(1) << 0, 2, 0;
    M1.row(2) << 1, 0, 2;
    M2.row(0) << 0, 1, 0;
    M2.row(1) << 1, 0, 1;
    M2.row(2) << 0, 1, 0;
    M.block<3, 3>(0, 0) = M1;
    M.block<3, 3>(3, 3) = M1;

    M.block<3, 3>(0, 3) = M2;
    M.block<3, 3>(3, 0) = M2;

    MatrixXf T;
    T = MatrixXf::Zero(3, 3);
    T.row(2) << 1, 1, 1;
    T.row(0) << coord(0, 0), coord(1, 0), coord(2, 0);
    T.row(1) << coord(0, 1), coord(1, 1), coord(2, 1);

    return D * N.transpose() * M * N * D / (24 * T.determinant());
};

inline MatrixXf MHFEM::Calculate_velocity_normal_vec_sep_edges(DFN::Domain dom,
                                                               size_t Frac_ID,
                                                               MatrixXf coord)
{

    MatrixXf Edge_outer_normal_2D = MatrixXf::Zero(3, 2);
    MatrixXf Edge_outer_normal_3D = MatrixXf::Zero(3, 3);

    RowVector2f node1 = coord.row(0),
                node2 = coord.row(1),
                node3 = coord.row(2);

    Edge_outer_normal_2D.row(0) << (node2 - node1)[1], -(node2 - node1)[0];
    Edge_outer_normal_2D.row(1) << (node3 - node2)[1], -(node3 - node2)[0];
    Edge_outer_normal_2D.row(2) << (node1 - node3)[1], -(node1 - node3)[0];

    for (size_t i = 0; i < 3; ++i)
    {
        std::vector<Vector3d> v_2d(1), v_3d(1);
        v_2d[0](0) = Edge_outer_normal_2D(i, 0);
        v_2d[0](1) = Edge_outer_normal_2D(i, 1);
        v_2d[0](2) = 0;

        DFN::Polygon_convex_3D poly{dom.Fractures[Frac_ID].Verts_trim};
        DFN::Rotate_to_horizontal R1{poly};
        R1.Rotate_back_without_z(v_2d, v_3d);

        Vector3f D;
        D << v_3d[0](0),
            v_3d[0](1),
            v_3d[0](2);

        D = D / (D.norm());

        Edge_outer_normal_3D.row(i) << D[0], D[1], D[2];
    }

    return Edge_outer_normal_3D;
};

void MHFEM::Matlab_plot(string FileKey_mat,
                        string FileKey_m,
                        DFN::Mesh_DFN_linear mesh,
                        DFN::Domain dom,
                        string same_dir)
{
    const char *filename = FileKey_mat.c_str();

    DFN::MATLAB_DATA_API M1_;
    M1_.Write_mat(filename, "w", 1, 1, 1, {0}, "nothing_");

    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i)
    {
        //cout << i << endl;
        MatrixXs ele_2D_frac = mesh.element_2D[i];

        vector<double> pData1(ele_2D_frac.rows() * 3);

        for (size_t j = 0; j < (size_t)ele_2D_frac.rows() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / ele_2D_frac.rows()); // column
            l = j % ele_2D_frac.rows();       // row

            pData1[j] = ele_2D_frac(l, k);
        }
        M1_.Write_mat(filename, "u", ele_2D_frac.rows() * 3,
                      ele_2D_frac.rows(), 3, pData1,
                      "element_2D_Frac_" + to_string(i + 1));

        //----------------------------------
        vector<double> pData2(mesh.coordinate_2D[i].rows() * 2);

        for (size_t j = 0; j < (size_t)mesh.coordinate_2D[i].rows() * 2; ++j)
        {
            size_t k, l;
            k = ceil(j / mesh.coordinate_2D[i].rows()); // column
            l = j % mesh.coordinate_2D[i].rows();       // row

            pData2[j] = mesh.coordinate_2D[i].coeffRef(l, k);
        }
        M1_.Write_mat(filename, "u", mesh.coordinate_2D[i].rows() * 2,
                      mesh.coordinate_2D[i].rows(), 2, pData2,
                      "coordinate_2D_Frac_" + to_string(i + 1));

        //----------------------------------
    }

    vector<double> pData4(mesh.coordinate_3D.rows() * 3);
    vector<double> pData5(mesh.element_3D.rows() * 3);
    vector<double> pData6(this->normal_vec_sep_edges.rows() * 3);
    vector<double> pData7(this->velocity_normal_scalar_sep_edges.rows());
    vector<double> pData8(this->pressure_interior_edge.rows());
    vector<double> pData9(this->pressure_eles.rows());

    for (size_t j = 0; j < (size_t)mesh.coordinate_3D.rows() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / mesh.coordinate_3D.rows()); // column
        l = j % mesh.coordinate_3D.rows();       // row

        pData4[j] = mesh.coordinate_3D(l, k);
        //cout << pData4[j] << endl;
    }

    for (size_t j = 0; j < (size_t)mesh.element_3D.rows() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / mesh.element_3D.rows()); // column
        l = j % mesh.element_3D.rows();       // row

        pData5[j] = mesh.element_3D(l, k);
    }

    for (size_t j = 0; j < (size_t)this->normal_vec_sep_edges.rows() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / this->normal_vec_sep_edges.rows()); // column
        l = j % this->normal_vec_sep_edges.rows();       // row

        pData6[j] = this->normal_vec_sep_edges(l, k);
    }

    for (size_t j = 0; j < (size_t)this->velocity_normal_scalar_sep_edges.rows(); ++j)
    {
        size_t k, l;
        k = ceil(j / this->velocity_normal_scalar_sep_edges.rows()); // column
        l = j % this->velocity_normal_scalar_sep_edges.rows();       // row

        pData7[j] = this->velocity_normal_scalar_sep_edges(l, k);
    }

    for (size_t j = 0; j < (size_t)pressure_interior_edge.rows(); ++j)
    {
        size_t k, l;
        k = ceil(j / pressure_interior_edge.rows()); // column
        l = j % pressure_interior_edge.rows();       // row

        pData8[j] = pressure_interior_edge(l, k);
        //cout << pData4[j] << endl;
    }

    for (size_t j = 0; j < (size_t)pressure_eles.rows(); ++j)
    {
        size_t k, l;
        k = ceil(j / pressure_eles.rows()); // column
        l = j % pressure_eles.rows();       // row

        pData9[j] = pressure_eles(l, k);
        //cout << pData4[j] << endl;
    }

    string FracJXY3D_s = "coordinate_3D";
    string FracJM_s = "element_3D";
    string ourter_normal = "normal_vec_sep_edges";
    string normal_velo = "velocity_normal_scalar_sep_edges";

    M1_.Write_mat(filename, "u", mesh.coordinate_3D.rows() * 3,
                  mesh.coordinate_3D.rows(), 3, pData4, FracJXY3D_s);
    M1_.Write_mat(filename, "u", mesh.element_3D.rows() * 3,
                  mesh.element_3D.rows(), 3, pData5, FracJM_s);
    M1_.Write_mat(filename, "u", this->normal_vec_sep_edges.rows() * 3,
                  this->normal_vec_sep_edges.rows(), 3, pData6, ourter_normal);
    M1_.Write_mat(filename, "u", this->velocity_normal_scalar_sep_edges.rows() * 1,
                  this->velocity_normal_scalar_sep_edges.rows(), 1, pData7, normal_velo);

    M1_.Write_mat(filename, "u", this->pressure_interior_edge.rows() * 1,
                  this->pressure_interior_edge.rows(), 1, pData8,
                  "pressure_interior_edge");

    M1_.Write_mat(filename, "u", this->pressure_eles.rows() * 1,
                  this->pressure_eles.rows(), 1, pData9,
                  "pressure_eles");

    //-------------------------
    vector<double> pData11(mesh.NUM_interior_edges * 2, 0);

    for (size_t i = 0; i < mesh.element_2D.size(); ++i)
        for (size_t j = 0; j < (size_t)mesh.element_2D[i].rows(); ++j)
            for (size_t k = 0; k < 3; ++k)
            {
                if (mesh.Interior_edgeNO[i](j, k).first == "interior")
                {
                    size_t globalID = mesh.Interior_edgeNO[i](j, k).second[0];

                    if (pData11[globalID - 1] == 0)
                    {
                        size_t node1 = mesh.element_2D[i](j, k);
                        size_t node2 = mesh.element_2D[i](j, (k + 1) % 3);

                        pData11[globalID - 1] = node1;
                        pData11[(globalID - 1) + mesh.NUM_interior_edges] = node2;
                    }
                }
            }
    M1_.Write_mat(filename, "u", mesh.NUM_interior_edges * 2,
                  mesh.NUM_interior_edges, 2, pData11, "Interior_edgeNO");

    // m file
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;";
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

    //oss << "figure(1); hold on\n";
    //oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', zeros(size(coordinate_3D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);\n";
    //oss << "view(3); hold on\n\n";

    oss << "element_2D_overall = [\n";
    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i)
        if (i != mesh.Frac_Tag.size() - 1)
            oss << "element_2D_Frac_" << to_string(i + 1) << endl;
        else
            oss << "element_2D_Frac_" << to_string(i + 1) << "];\n";
    oss << "pressure_x = []; pressure_y = []; pressure_z = []; pressure_ = pressure_interior_edge; tmp = 1;\n";

    oss << "for i = 1:size(Interior_edgeNO, 1)\n";
    oss << "\tnode1 = Interior_edgeNO(i, 1); node2 = Interior_edgeNO(i, 2);\n";
    oss << "\tcoord_t = 0.5 * (coordinate_3D(node1, :) + coordinate_3D(node2, :));\n";
    oss << "\tpressure_x(i, 1) = coord_t(1, 1);\n";
    oss << "\tpressure_y(i, 1) = coord_t(1, 2);\n";
    oss << "\tpressure_z(i, 1) = coord_t(1, 3);\n";
    oss << "end\n";

    oss << "m_t = size(pressure_x, 1);\ncenter_ele = [];\n";

    oss << "for i = 1:size(element_2D_overall, 1)\n";
    oss << "\tnode_t = element_2D_overall(i, :);\n";
    oss << "\tcenter_ele(i, :) = (coordinate_3D(node_t(1), :) + coordinate_3D(node_t(2), :) + coordinate_3D(node_t(3), :)) * (1/3);\n";
    oss << "end\n";

    oss << "pressure_x = [pressure_x; center_ele(:, 1)];\n";
    oss << "pressure_y = [pressure_y; center_ele(:, 2)];\n";
    oss << "pressure_z = [pressure_z; center_ele(:, 3)];\n";
    oss << "pressure_ = [pressure_; pressure_eles];\n\n";

    if (same_dir == "NO")
        oss << "pressure_vert = griddata(pressure_x, pressure_y, pressure_z, pressure_, coordinate_3D(:, 1), coordinate_3D(:, 2), coordinate_3D(:, 3), 'nearest');\n";
    else if (same_dir == "x")
        oss << "pressure_vert = griddata(pressure_y, pressure_z, pressure_, coordinate_3D(:, 2), coordinate_3D(:, 3), 'nearest');\n";
    else if (same_dir == "y")
        oss << "pressure_vert = griddata(pressure_x, pressure_z, pressure_, coordinate_3D(:, 1), coordinate_3D(:, 3), 'nearest');\n";
    else if (same_dir == "z")
        oss << "pressure_vert = griddata(pressure_x, pressure_y, pressure_, coordinate_3D(:, 1), coordinate_3D(:, 2), 'nearest');\n";

    oss << "figure(1); hold on\n";
    oss << "patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_vert, 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); colorbar; view(3)\n";
    oss << "caxis([" << this->outlet_p << ", " << this->inlet_p << "]);\n";
    //------------------------
    oss << "\n\ncenter = zeros(size(element_2D_overall, 1)*3, 3); edgeno = 1;\n";
    oss << "for i = 1:size(element_2D_overall, 1)\n";
    oss << "\tfor j = 1:3\n";
    oss << "\t\tnode1 = element_2D_overall(i,j); node2 = element_2D_overall(i, mod(j, 3) + 1);\n";
    oss << "\t\tcenter(edgeno, :) = 0.5 * (coordinate_3D(node1, :) + coordinate_3D(node2, :)); edgeno = edgeno + 1;\n";
    oss << "\tend\n";
    oss << "end\n";
    oss << "normal_vec_sep_edges = normal_vec_sep_edges.* velocity_normal_scalar_sep_edges;\n";
    oss << "hold on;\nquiver3(center(:, 1), center(:, 2), center(:, 3), ...\n";
    oss << "\tnormal_vec_sep_edges(:, 1), normal_vec_sep_edges(:, 2), normal_vec_sep_edges(:, 3), 3, 'linewidth', 1.3, 'color', 'r');\n";

    oss.close();
};

inline double MHFEM::The_inlet_of_all_elements(DFN::Mesh_DFN_linear mesh)
{
    size_t _nodeid_ = 0;
    double inlet_all_ele = 0;

    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i) // Frac
    {
        for (size_t j = 0; j < (size_t)mesh.element_2D[i].rows(); ++j) // elements
        {
            for (size_t k = 0; k < 3; ++k)
            {
                if (velocity_normal_scalar_sep_edges(_nodeid_, 0) < 0)
                {
                    size_t node_a = mesh.element_2D[i](j, k),
                           node_b = mesh.element_2D[i](j, (k + 1) % 3);

                    double len = (mesh.coordinate_3D.row(node_a - 1) - mesh.coordinate_3D.row(node_b - 1)).norm();

                    inlet_all_ele += len * velocity_normal_scalar_sep_edges(_nodeid_, 0);
                }
                _nodeid_++;
            }
        }
    }
    return inlet_all_ele;
};

inline size_t MHFEM::Element_ID(size_t i, size_t j, DFN::Mesh_DFN_linear mesh)
{
    size_t ID = 0;

    if (i == 0)
        return ID + j;

    for (size_t y = 0; y < i; ++y)
        ID += mesh.element_2D[y].rows();

    return ID + j;
};
}; // namespace DFN