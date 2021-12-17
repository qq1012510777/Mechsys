#pragma once
#include "../Error_throw/Error_throw.h"
#include "../Geometry_H/Normal_vector_2D.h"
#include "../MATLAB_DATA_API/MATLAB_DATA_API.h"

#include "../Mesh_H/Mesh_DFN_linear.h"
#include "../Using_UMFPACK/Using_UMFPACK.h"

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
    size_t dir = 0;

    MatrixXf velocity_normal_scalar_sep_edges;
    MatrixXf pressure_interior_edge;
    MatrixXf pressure_eles;

    MatrixXf normal_vec_sep_edges;

public:
    MHFEM(DFN::Mesh_DFN_linear mesh,
          DFN::Domain dom,
          double inlet_p_,
          double outlet_p_,
          size_t i /*direction*/);

    void Matlab_plot(string FileKey_mat,
                     string FileKey_m,
                     DFN::Mesh_DFN_linear mesh,
                     DFN::Domain dom);

private:
    void Assemble_and_solve(DFN::Mesh_DFN_linear mesh,
                            DFN::Domain dom);

    MatrixXf Calculate_velocity_normal_vec_sep_edges(DFN::Domain dom,
                                                     size_t Frac_ID,
                                                     MatrixXf coord);

private:
    MatrixXf stimaB(MatrixXf coord);
};

inline MHFEM::MHFEM(DFN::Mesh_DFN_linear mesh,
                    DFN::Domain dom,
                    double inlet_p_,
                    double outlet_p_,
                    size_t i /*direction*/)
{
    dir = i;
    inlet_p = inlet_p_;
    outlet_p = outlet_p_;
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
    size_t NUM_sep_edges = mesh.element_3D.rows() * 3,
           NUM_eles = mesh.element_3D.rows(),
           NUM_glob_interior_edges = mesh.NUM_interior_edges;

    size_t Dim = NUM_sep_edges + NUM_eles + NUM_glob_interior_edges;

    normal_vec_sep_edges.resize(NUM_sep_edges, 3);

    MatrixXf A = MatrixXf::Zero(NUM_sep_edges, NUM_sep_edges);
    MatrixXf B = MatrixXf::Zero(NUM_sep_edges, NUM_eles);
    MatrixXf C = MatrixXf::Zero(NUM_glob_interior_edges, NUM_sep_edges);

    size_t tmp_eleNO = 0;

    vector<pair<size_t, double>> Inlet_edge_sep_NO,
        Outlet_edge_sep_NO;
    vector<size_t> Neumann_edge_sep_NO;

    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i) // Frac
    {
        size_t Frac_Tag = mesh.Frac_Tag[i];

        double Kper = dom.Fractures[Frac_Tag].Conductivity;
        for (size_t j = 0; j < (size_t)mesh.element_2D[i].rows(); ++j) // elements
        {

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

            MatrixXf outnormal_vec = Calculate_velocity_normal_vec_sep_edges(dom, Frac_Tag, coord);
            normal_vec_sep_edges.block(edge1 - 1, 0, 3, 3) = outnormal_vec;

            //vector<size_t> I{edge2 - 1, edge3 - 1, edge1 - 1};
            vector<size_t> I{edge2 - 1, edge3 - 1, edge1 - 1};

            for (size_t ik = 0; ik < 3; ++ik)
                for (size_t jk = 0; jk < 3; ++jk)
                {
                    A(I[ik], I[jk]) += round(((-1.0 / Kper) * A_(ik, jk)), 4);
                }
            /*
            if(node1 == 9 && node2 == 4 && node3 == 2)
            {
                cout << "----------\n";
                cout << A_ << endl;
                cout << "----------\n";
            }
            */

            Vector3f B_;
            B_ << round((coord.row(2) - coord.row(1)).norm(), 4),
                round((coord.row(0) - coord.row(2)).norm(), 4),
                round((coord.row(1) - coord.row(0)).norm(), 4);

            for (size_t ik = 0; ik < 3; ++ik)
                B(I[ik], tmp_eleNO) = round(B_(ik, 0), 4);

            tmp_eleNO++;

            // how to know if an edge is an interior edge
            //cout << 1 << endl;
            for (size_t ik = 0; ik < 3; ++ik)
            {
                size_t ik_r = (ik + 1) % 3; // local edgeNO

                size_t node_a = mesh.element_2D[i](j, ik_r);
                size_t node_b = mesh.element_2D[i](j, (ik_r + 1) % 3);

                size_t node_A = node_a < node_b ? node_a : node_b,
                       node_B = node_a > node_b ? node_a : node_b;

                std::map<pair<size_t, size_t>, size_t>::iterator ity;
                ity = mesh.Interior_edgeNO.find(make_pair(node_A, node_B));

                double len = (mesh.coordinate_3D.row(node_A - 1) - mesh.coordinate_3D.row(node_B - 1)).norm();

                if (ity != mesh.Interior_edgeNO.end())
                {
                    //this edge is an interior edge
                    size_t global_inter_edge_NO = ity->second;
                    C(global_inter_edge_NO - 1, I[ik]) -= round(len, 4);
                    continue;
                }

                //------------------------
                //Dirichlet inlet
                std::set<pair<size_t, size_t>>::iterator it_inlet;
                it_inlet = mesh.Inlet_edges.find(make_pair(node_A, node_B));
                if (it_inlet != mesh.Inlet_edges.end())
                {
                    Inlet_edge_sep_NO.push_back(std::make_pair(I[ik] + 1, len));
                    continue;
                }

                //Dirichlet outlet
                std::set<pair<size_t, size_t>>::iterator it_outlet;
                it_outlet = mesh.Outlet_edges.find(make_pair(node_A, node_B));
                if (it_outlet != mesh.Outlet_edges.end())
                {
                    Outlet_edge_sep_NO.push_back(std::make_pair(I[ik] + 1, len));
                    continue;
                }

                //Neumann edges
                std::set<pair<size_t, size_t>>::iterator Neumann_edges_;
                Neumann_edges_ = mesh.Neumann_edges.find(make_pair(node_A, node_B));
                if (Neumann_edges_ != mesh.Neumann_edges.end())
                {
                    Neumann_edge_sep_NO.push_back(I[ik] + 1);
                    continue;
                }

                throw Error_throw_pause("An edge must has an attribute\n");
            }
        }
    }

    MatrixXf K = MatrixXf::Zero(Dim, Dim);
    MatrixXf b = MatrixXf::Zero(Dim, 1);

    K.block(0, 0, NUM_sep_edges, NUM_sep_edges) = A;
    K.block(0, NUM_sep_edges, NUM_sep_edges, NUM_eles) = B;
    K.block(0, NUM_sep_edges + NUM_eles, NUM_sep_edges, NUM_glob_interior_edges) = C.transpose();
    K.block(NUM_sep_edges, 0, NUM_eles, NUM_sep_edges) = B.transpose();
    K.block(NUM_sep_edges + NUM_eles, 0, NUM_glob_interior_edges, NUM_sep_edges) = C;

    // Dirichlet boundary condition
    for (size_t i = 0; i < Inlet_edge_sep_NO.size(); ++i)
    {
        size_t SEP_edgeNO = Inlet_edge_sep_NO[i].first;
        double len = Inlet_edge_sep_NO[i].second;

        double P_D = this->inlet_p;
        b(SEP_edgeNO - 1, 0) += P_D * len;
    }
    for (size_t i = 0; i < Outlet_edge_sep_NO.size(); ++i)
    {
        size_t SEP_edgeNO = Outlet_edge_sep_NO[i].first;
        double len = Outlet_edge_sep_NO[i].second;

        double P_D = this->outlet_p;
        b(SEP_edgeNO - 1, 0) += P_D * len;
    }

    //--------Neumann
    for (size_t i = 0; i < Neumann_edge_sep_NO.size(); ++i)
    {
        size_t SEP_edgeNO = Neumann_edge_sep_NO[i];
        K.col(SEP_edgeNO - 1) = VectorXf::Zero(Dim, 1);
        K.row(SEP_edgeNO - 1) = RowVectorXf::Zero(1, Dim);
        K(SEP_edgeNO - 1, SEP_edgeNO - 1) = 1;
        b(SEP_edgeNO - 1, 0) = 0;
    }

    //MatrixXf x = K.inverse() * b;
    //cout << K << endl;
    //cout << "\nb\n";
    //cout << b << endl;
    //cout << "\nx\n";
    //cout << x << endl;

    B.resize(NUM_eles, NUM_sep_edges);

    A = K.block(0, 0, NUM_sep_edges, NUM_sep_edges);
    B = K.block(NUM_sep_edges, 0, NUM_eles, NUM_sep_edges);
    C = K.block(NUM_sep_edges + NUM_eles, 0, NUM_glob_interior_edges, NUM_sep_edges);

    MatrixXf g = b.block(0, 0, NUM_sep_edges, 1);
    MatrixXf f = b.block(NUM_sep_edges, 0, NUM_eles, 1);

    MatrixXf A_inv = A.inverse(); //------

    MatrixXf C_tps = C.transpose();
    MatrixXf B_tps = B.transpose();
    MatrixXf Wq = B * A_inv;
    MatrixXf U = (Wq * B_tps).inverse();
    MatrixXf Re = C * A_inv;
    MatrixXf Eq = Re * B_tps;
    MatrixXf Sd = Wq * C_tps;

    MatrixXf D = Re * C_tps - Eq * U * Sd;
    MatrixXf r = Re * g + Eq * U * (f - Wq * g);

    SparseMatrix<float> D_ = D.sparseView();
    ConjugateGradient<SparseMatrix<float>, Lower | Upper> cg;
    cg.compute(D_);
    pressure_interior_edge = cg.solve(r);

    pressure_eles = U * (Wq * g - Sd * pressure_interior_edge - f);
    velocity_normal_scalar_sep_edges = A_inv * (g - B_tps * pressure_eles - C_tps * pressure_interior_edge);

    for (size_t i = 0; i < Inlet_edge_sep_NO.size(); ++i)
    {
        size_t sep_EDGE_no = Inlet_edge_sep_NO[i].first;
        double len = Inlet_edge_sep_NO[i].second;

        double veloc_length = abs(velocity_normal_scalar_sep_edges(sep_EDGE_no - 1, 0) * len);
        this->Q_in += veloc_length;
    }
    for (size_t i = 0; i < Outlet_edge_sep_NO.size(); ++i)
    {
        size_t sep_EDGE_no = Outlet_edge_sep_NO[i].first;
        double len = Outlet_edge_sep_NO[i].second;

        double veloc_length = abs(velocity_normal_scalar_sep_edges(sep_EDGE_no - 1, 0) * len);
        this->Q_out += veloc_length;
    }
};

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
                        DFN::Domain dom)
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
    vector<double> pData11(mesh.Interior_edgeNO.size() * 2);

    for (std::map<pair<size_t, size_t>, size_t>::iterator its = mesh.Interior_edgeNO.begin();
         its != mesh.Interior_edgeNO.end(); its++)
    {
        size_t yw = its->second;
        size_t node1 = its->first.first;
        pData11[yw - 1] = node1;
    }

    for (std::map<pair<size_t, size_t>, size_t>::iterator its = mesh.Interior_edgeNO.begin();
         its != mesh.Interior_edgeNO.end(); its++)
    {
        size_t yw = its->second;
        size_t node1 = its->first.second;
        pData11[yw - 1 + mesh.Interior_edgeNO.size()] = node1;
    }

    M1_.Write_mat(filename, "u", mesh.Interior_edgeNO.size() * 2,
                  mesh.Interior_edgeNO.size(), 2, pData11, "Interior_edgeNO");

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

    oss << "pressure_vert = griddata(pressure_x, pressure_y, pressure_z, pressure_, coordinate_3D(:, 1), coordinate_3D(:, 2), coordinate_3D(:, 3), 'nearest');\n";
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

}; // namespace DFN