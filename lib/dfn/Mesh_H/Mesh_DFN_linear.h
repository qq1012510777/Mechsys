#pragma once
#include "../DFN_H/Domain_WL.h"
#include "../Eigen_API/Eigen_API.h"
#include "../Geometry_H/If_skinny_triangle.h"
#include "../Geometry_H/Point_3D.h"
#include "../Geometry_H/Triangle_orientation.h"
#include "../Remove_fractures/Remove_fractures.h"
#include "Eigen/Sparse"
#include "mat.h"
#include <map>

//
#include <gmsh.h>

namespace DFN
{
class Mesh_DFN_linear
{
public:
    bool mesh_state = true;
    std::vector<size_t> Frac_Tag;

    std::vector<Sp_d_mat> coordinate_2D;
    std::vector<MatrixXs> element_2D;

    MatrixXd coordinate_3D;
    vector<Vector7b> _Pnt_attri;
    //bool If_model_right = false;
    //bool If_model_left = false;
    //bool If_model_back = false;
    //bool If_model_front = false;
    //bool If_model_bottom = false;
    //bool If_model_top = false;
    //bool If_frac_bound = false;

    MatrixXs element_3D;
    //vector<Sp_s_mat> nodes2edge;
    //vector<Sp_s_mat> nodes2ele;
    //vector<size_t> Num_edges_frac;
    //vector<MatrixXs> edge2element;
    //vector<MatrixXs> Neumann_edges;
    //vector<MatrixXs> Dirichlet_edges;
    //vector<MatrixXs> Dirichlet_attri;
    //vector<MatrixXs> Interior_edges;

    //vector<MatrixX_st_v2d> Interior_edgeNO;
    //vector<vector<Vector3d>> Interior_edgeNO_to_ele_edge;
    //std::map<pair<size_t, size_t>, size_t> Inlet_edges;
    //std::map<pair<size_t, size_t>, size_t> Outlet_edges;
    //std::map<pair<size_t, size_t>, size_t> Neumann_edges;
    //SparseMatrix<size_t> Interior_edgeNO;

    vector<MatrixX_st_v2d> Edge_attri_and_NO;

    size_t NUM_interior_edges = 0;
    size_t NUM_inlet_edges = 0;
    size_t NUM_outlet_edges = 0;
    size_t NUM_neumann_edges = 0;
    size_t dir_ = 2;

    SparseMatrix<double> Edge_length;
    double mean_edge_length = 0;
    double max_edge_length = 0;
    double min_edge_length = 1e5;

    size_t Nproc_global = 1;

public:
    Mesh_DFN_linear();
    Mesh_DFN_linear(DFN::Domain dom,
                    const double min_ele_edge,
                    const double max_ele_edge,
                    size_t dir,
                    size_t Nproc);

    void Matlab_plot(string FileKey_mat,
                     string FileKey_m,
                     DFN::Domain dom);

private:
    void Generate_2D_nodes(DFN::Domain dom);
    void Identify_point_attribute(DFN::Domain dom);
    void Numbering_edges(DFN::Domain dom);
    void Modify_the_triangle_orientation();
    void Identify_ele_fracNO(DFN::Domain dom);

private:
    bool If_interior_edge(size_t node1, size_t node2);
    bool If_two_pnts_Dirchlet(size_t node1, size_t node2, string &A);
    bool If_two_pnts_Neumann(size_t node1, size_t node2);
    size_t Element_ID(size_t i, size_t j);
    SparseMatrix<double> Node2element(size_t i);
};

inline Mesh_DFN_linear::Mesh_DFN_linear(){};

inline Mesh_DFN_linear::Mesh_DFN_linear(DFN::Domain dom,
                                        const double min_ele_edge,
                                        const double max_ele_edge,
                                        size_t dir,
                                        size_t Nproc)
{
    this->dir_ = dir;
    this->Nproc_global = Nproc;
    //---------------------------------------
    if (dom.Percolation_cluster[dir].size() > 1)
    {
        string AS = "The DFN has more than one percolating cluster\n";
        //throw Error_throw_ignore(AS);
    }

    if (dom.Percolation_cluster[dir].size() == 0)
    {
        string AS = "The DFN has no one percolating cluster\n";
        throw Error_throw_ignore(AS);
    }

    try
    {
        size_t VertsPntID = 1;
        size_t LineID = 1;
        size_t CurveLoopID = 1;
        size_t SurfaceID = 1;

        //cout << "mesh" << endl;
        gmsh::initialize();
        //cout << "mesh init" << endl;
        //gmsh::option::setNumber("General.NumThreads", Nproc);
        gmsh::option::setNumber("General.Verbosity", 2); // default level is 5
        gmsh::model::add("t2");
        gmsh::model::mesh::clear();
        //cout << "size: " << dom.Percolation_cluster[dir].size() << endl;
        //cout << "status: " << dom.Percolation_status[dir] << endl;

        size_t iF_only_one_frac = false;
        for (size_t i = 0; i < dom.Percolation_cluster[dir].size(); ++i)
        {
            size_t ClusterID = dom.Percolation_cluster[dir][i];

            if (dom.Listofclusters[ClusterID].size() == 1 && dom.Percolation_cluster[dir].size() == 1)
                iF_only_one_frac = true;

            DFN::Remove_fractures Remo;

            vector<size_t> LISTCOLUSTER = dom.Listofclusters[ClusterID];

            /*
            for (size_t j = 0; j < LISTCOLUSTER.size(); ++j)
                cout << "[" << LISTCOLUSTER[j] << ", " << j + 1 << "], ";
            cout << endl;
            */

            for (size_t j = 0;; ++j)
            {
                bool ses1 = Remo.Remove_dead_fractures(LISTCOLUSTER, dom, this->dir_);
                if (ses1 == false)
                    break;
            }

            /*
            for (size_t j = 0; j < LISTCOLUSTER.size(); ++j)
                cout << LISTCOLUSTER[j] << ", ";
            cout << endl;
            */

            for (size_t j = 0; j < LISTCOLUSTER.size(); ++j)
            {
                size_t FracID = LISTCOLUSTER[j];

                this->Frac_Tag.push_back(FracID);

                std::vector<int> Pointloop(dom.Fractures[FracID].Verts_trim.size());

                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                    Pointloop[k] = VertsPntID + k;

                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                {

                    gmsh::model::occ::addPoint(dom.Fractures[FracID].Verts_trim[k](0),
                                               dom.Fractures[FracID].Verts_trim[k](1),
                                               dom.Fractures[FracID].Verts_trim[k](2),
                                               0,
                                               VertsPntID);
                    VertsPntID++;
                }

                std::vector<int> curveloop(dom.Fractures[FracID].Verts_trim.size());

                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                    curveloop[k] = LineID + k;

                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                {
                    gmsh::model::occ::addLine(Pointloop[k],
                                              Pointloop[(k + 1) % dom.Fractures[FracID].Verts_trim.size()],
                                              LineID);
                    LineID++;
                }
                gmsh::model::occ::addCurveLoop(curveloop, CurveLoopID);

                std::vector<int> surfaceloop = {(int)CurveLoopID};
                gmsh::model::occ::addPlaneSurface(surfaceloop, (int)SurfaceID);

                CurveLoopID++;
                SurfaceID++;
            }
        }
        gmsh::model::occ::synchronize();

        if (iF_only_one_frac == false)
        {

            std::vector<std::pair<int, int>> input_entity(SurfaceID - 1);

            for (size_t i = 0; i < input_entity.size(); ++i)
                input_entity[i] = std::make_pair(2, i + 1);

            std::vector<std::pair<int, int>> out;
            std::vector<std::vector<std::pair<int, int>>> outmap;

            gmsh::model::occ::fragment(input_entity, input_entity, out, outmap);
            gmsh::model::occ::synchronize();
        }
        gmsh::option::setNumber("Mesh.MeshSizeMin", min_ele_edge);
        gmsh::option::setNumber("Mesh.MeshSizeMax", max_ele_edge);

        gmsh::option::setNumber("Mesh.Algorithm", 5);

        //std::cout << "\033[31mstart meshing;\n\033[0m";
        gmsh::model::mesh::generate(2);
        //std::cout << "\033[31mfinish meshing;\n\033[0m";
        //--------

        //-----------
        std::vector<std::size_t> nodes;
        std::vector<double> coord, coordParam;
        gmsh::model::mesh::getNodes(nodes, coord, coordParam);
        size_t NUM_nodes = coord.size() / 3;

        coordinate_3D.resize(NUM_nodes, 3);

        for (size_t i = 0; i < coord.size(); i += 3)
            coordinate_3D.row(i / 3) << coord[i], coord[i + 1], coord[i + 2];

        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 2, -1);
        size_t NUM_ele = elemNodeTags[0].size() / 3;
        this->element_3D.resize(NUM_ele, 3);

        vector<int> row_to_revomve;
        for (size_t i = 0; i < elemNodeTags[0].size(); i += 3)
        {
            //cout << (size_t)elemNodeTags[0][i] << ", " << (size_t)elemNodeTags[0][i + 1] << ", " << (size_t)elemNodeTags[0][i + 2] << endl;
            this->element_3D.row(i / 3) << (size_t)elemNodeTags[0][i], (size_t)elemNodeTags[0][i + 1], (size_t)elemNodeTags[0][i + 2];

            vector<RowVector3d> coord_(3);

            coord_[0] = this->coordinate_3D.row((size_t)elemNodeTags[0][i] - 1);
            coord_[1] = this->coordinate_3D.row((size_t)elemNodeTags[0][i + 1] - 1);
            coord_[2] = this->coordinate_3D.row((size_t)elemNodeTags[0][i + 2] - 1);

            //cout << (size_t)elemNodeTags[0][i] << ", " << (size_t)elemNodeTags[0][i + 1] << ",  " << (size_t)elemNodeTags[0][i + 2] << "\n\t";
            DFN::If_skinny_triangle SKINNY{coord_};

            if (SKINNY.If_skinny == true)
            {
                row_to_revomve.push_back(i / 3);
                //cout << "element TAG: " << i / 3 + 1 << endl;
            }
        }

        //cout << 1 << endl;
        if (row_to_revomve.size() >= 1)
        {
            MatrixXs ATR = removeMatrixXsRow(this->element_3D, row_to_revomve);
            this->element_3D.resize(this->element_3D.rows() - row_to_revomve.size(), 3);
            this->element_3D = ATR;
        }
        //cout << 2 << endl;
        //------------2D

        //-----------
        //cout << "Identify_ele_fracNO\n";
        this->Identify_ele_fracNO(dom);
        //cout << "Generate_2D_nodes\n";
        this->Generate_2D_nodes(dom);
        //cout << "Modify_the_triangle_orientation\n";
        this->Modify_the_triangle_orientation();
        this->Identify_point_attribute(dom);
        //cout << this->element_2D[0] << endl;
        //cout << "-------\n";
        //cout << "Numbering_edges\n";
        this->Numbering_edges(dom);

        /*
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTags;

        gmsh::model::mesh::getElements(elementTypes,
                                       elementTags,
                                       nodeTags,
                                       2,
                                       4);
        cout << nodeTags[0].size() / 3 << endl;
        */

        //gmsh::fltk::run();
        gmsh::model::mesh::clear();
        gmsh::clear();
        gmsh::finalize();
        //cout << "mesh finished\n";
    }
    catch (...)
    {
        this->mesh_state = false;
        gmsh::clear();
        gmsh::finalize();
    }
};

inline void Mesh_DFN_linear::Generate_2D_nodes(DFN::Domain dom)
{
    this->coordinate_2D.resize(this->Frac_Tag.size());
    size_t Nproc_ = this->Nproc_global;

#pragma omp parallel for schedule(dynamic) num_threads(Nproc_)
    for (size_t i = 0; i < this->Frac_Tag.size(); ++i)
    {
        Eigen::SparseMatrix<double> m1(this->coordinate_3D.rows(), 2);
        m1.reserve(VectorXi::Constant(2, this->coordinate_3D.rows()));

        VectorXd PNT_INDICATOR = Eigen::VectorXd::Zero(this->coordinate_3D.rows());

        size_t Frac_Tag = this->Frac_Tag[i];
        Vector3d add_;
        add_ << random_double(1, 10), random_double(1, 10), random_double(1, 10);
        MatrixXd A = dom.Fractures[Frac_Tag].Rotation_Matrix_3DTo2D(add_);

        for (size_t j = 0; j < (size_t)element_2D[i].rows(); ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                //cout << 1 << endl;
                size_t PNT_ID = element_2D[i](j, k) - 1;
                //cout << "PNT_ID " << PNT_ID << endl;
                if (PNT_INDICATOR[PNT_ID] == 0)
                {
                    /////cout << 1.1 << endl;
                    PNT_INDICATOR[PNT_ID] = 1;
                    //cout << 1.2 << endl;

                    Vector4d PNT_RR, PNT_WW;
                    PNT_RR << (double)this->coordinate_3D(PNT_ID, 0) + add_[0],
                        (double)this->coordinate_3D(PNT_ID, 1) + add_[1],
                        (double)this->coordinate_3D(PNT_ID, 2) + add_[2],
                        0;

                    //cout << 1.3 << endl;
                    PNT_WW = A * PNT_RR;

                    //cout << 1.5 << endl;
                    if (abs(PNT_WW[2]) > 1e-3)
                    {
                        cout << PNT_WW.transpose() << endl;
                        cout << "Rotate pnts to 2D failed! In class 'Mesh_DFN_linear'!\n";
                        throw Error_throw_ignore("Rotate pnts to 2D failed! In class 'Mesh_DFN_linear'!\n");
                    }
                    m1.insert(PNT_ID, 0) = PNT_WW[0];
                    m1.insert(PNT_ID, 1) = PNT_WW[1];
                }
                //cout << 2 << endl;
            }
        }
        m1.makeCompressed();
        this->coordinate_2D[i] = m1;
        //cout << m1.rows() << ", " << m1.cols() << endl;
    }
};

inline void Mesh_DFN_linear::Modify_the_triangle_orientation()
{

    std::vector<bool> If_clockwise(this->Frac_Tag.size());

    size_t Nproc_ = this->Nproc_global;
#pragma omp parallel for schedule(dynamic) num_threads(Nproc_)
    for (size_t i = 0; i < If_clockwise.size(); ++i)
    {

        std::vector<Vector2d> Triangle_(3);
        size_t pntID = this->element_2D[i](0, 0) - 1;

        Triangle_[0] << this->coordinate_2D[i].coeffRef(pntID, 0),
            this->coordinate_2D[i].coeffRef(pntID, 1);

        pntID = this->element_2D[i](0, 1) - 1;

        Triangle_[1] << this->coordinate_2D[i].coeffRef(pntID, 0),
            this->coordinate_2D[i].coeffRef(pntID, 1);

        pntID = this->element_2D[i](0, 2) - 1;
        Triangle_[2] << this->coordinate_2D[i].coeffRef(pntID, 0),
            this->coordinate_2D[i].coeffRef(pntID, 1);

        DFN::Triangle_orientation ori{Triangle_};

        if (ori.If_clockwise == true)
        {
            for (size_t j = 0; j < (size_t)this->element_2D[i].rows(); ++j)
            {
                size_t node1_ = this->element_2D[i](j, 0),
                       node2_ = this->element_2D[i](j, 2),
                       node3_ = this->element_2D[i](j, 1);

                this->element_2D[i].row(j) << node1_,
                    node2_,
                    node3_;
            }
        }
    }
};

inline void Mesh_DFN_linear::Numbering_edges(DFN::Domain dom)
{
    Edge_attri_and_NO.resize(element_2D.size());

    SparseMatrix<size_t> Shared_edge_global_NO(this->coordinate_3D.rows(), this->coordinate_3D.rows());
    Shared_edge_global_NO.reserve(VectorXi::Constant(this->coordinate_3D.rows(), 4));

    Edge_length.resize(this->coordinate_3D.rows(), this->coordinate_3D.rows());
    Edge_length.reserve(VectorXi::Constant(this->coordinate_3D.rows(), 1));

    size_t edge_shared = 1;

    size_t Sep_edge_NO_in = 1;
    size_t Sep_edge_NO_out = 1;
    size_t Sep_edge_NO_neumann = 1;

    for (size_t i = 0; i < element_2D.size(); ++i)
    {
        SparseMatrix<double> node2element_local = Node2element(i);
        Edge_attri_and_NO[i].resize(element_2D[i].rows(), 3);
        for (size_t j = 0; j < (size_t)element_2D[i].rows(); ++j)
        {
            size_t tmp_ele_NO = Element_ID(i, j);

            for (size_t k = 0; k < 3; ++k)
            {
                size_t Sep_NO = tmp_ele_NO * 3 + k + 1;

                Edge_attri_and_NO[i](j, k).first = "nothing";
                Edge_attri_and_NO[i](j, k).second << 0, 0;

                size_t node1 = element_2D[i](j, k),
                       node2 = element_2D[i](j, (k + 1) % 3);

                if (node2element_local.coeffRef(node1 - 1, node2 - 1) != 0 && node2element_local.coeffRef(node2 - 1, node1 - 1) != 0)
                {
                    //shared edge
                    Edge_attri_and_NO[i](j, k).first = "interior";
                    Edge_attri_and_NO[i](j, k).second << 0, 0;

                    if (Shared_edge_global_NO.coeffRef(node1 - 1, node2 - 1) == 0)
                    {
                        Edge_attri_and_NO[i](j, k).second << edge_shared, Sep_NO;
                        Shared_edge_global_NO.coeffRef(node1 - 1, node2 - 1) = edge_shared;
                        Shared_edge_global_NO.coeffRef(node2 - 1, node1 - 1) = edge_shared;
                        edge_shared++;
                    }
                    else
                        Edge_attri_and_NO[i](j, k).second << Shared_edge_global_NO.coeffRef(node1 - 1, node2 - 1), Sep_NO;
                }
                else
                {
                    // boundary edge
                    string A_t;
                    bool ty_1 = If_two_pnts_Dirchlet(node1, node2, A_t);

                    if (ty_1 == true)
                    {
                        if (A_t == "in")
                        {
                            Edge_attri_and_NO[i](j, k).first = "in";
                            Edge_attri_and_NO[i](j, k).second << Sep_NO, Sep_edge_NO_in;
                            Sep_edge_NO_in++;
                        }
                        else
                        {
                            Edge_attri_and_NO[i](j, k).first = "out";
                            Edge_attri_and_NO[i](j, k).second << Sep_NO, Sep_edge_NO_out;
                            Sep_edge_NO_out++;
                        }
                    }
                    else
                    {
                        Edge_attri_and_NO[i](j, k).first = "neumann";
                        Edge_attri_and_NO[i](j, k).second << Sep_NO, Sep_edge_NO_neumann;
                        Sep_edge_NO_neumann++;
                    }
                }

                double len = (this->coordinate_3D.row(node1 - 1) -
                              this->coordinate_3D.row(node2 - 1))
                                 .norm();

                std::pair<size_t, size_t> At_ = std::make_pair(node1 < node2 ? node1 : node2,
                                                               node1 > node2 ? node1 : node2);

                this->Edge_length.coeffRef(At_.first - 1, At_.second - 1) = len;
            }
        }
    }

    this->Edge_length.makeCompressed();

    NUM_interior_edges = edge_shared - 1;
    NUM_inlet_edges = Sep_edge_NO_in - 1;
    NUM_outlet_edges = Sep_edge_NO_out - 1;
    NUM_neumann_edges = Sep_edge_NO_neumann - 1;

    double len_totoal = 0;
    for (int k = 0; k < Edge_length.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(Edge_length, k); it; ++it)
        {
            double len = Edge_length.coeffRef(it.row(), it.col());

            this->max_edge_length = max_edge_length > len ? max_edge_length : len;
            this->min_edge_length = min_edge_length < len ? min_edge_length : len;
            len_totoal += len;
        }
    }

    this->mean_edge_length = len_totoal / (NUM_interior_edges + NUM_inlet_edges + NUM_outlet_edges + NUM_neumann_edges);
}

inline bool Mesh_DFN_linear::If_two_pnts_Neumann(size_t node1, size_t node2)
{
    size_t a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;
    if (this->dir_ == 0)
        a1 = 2, a2 = 3, a3 = 4, a4 = 5, a5 = 6;
    else if (this->dir_ == 1)
        a1 = 0, a2 = 1, a3 = 4, a4 = 5, a5 = 6;
    else if (this->dir_ == 2)
        a1 = 0, a2 = 1, a3 = 2, a4 = 3, a5 = 6;

    vector<size_t> A = {a1, a2, a3, a4, a5};

    bool Neumann1 = false, Neumann2 = false;

    for (size_t i = 0; i < 5; ++i)
        if (_Pnt_attri[node1 - 1][A[i]] == true)
        {
            Neumann1 = true;
            break;
        }

    for (size_t i = 0; i < 5; ++i)
        if (_Pnt_attri[node2 - 1][A[i]] == true)
        {
            Neumann2 = true;
            break;
        }

    if (Neumann1 == true && Neumann2 == true)
        return true;
    else
        return false;
};

inline bool Mesh_DFN_linear::If_two_pnts_Dirchlet(size_t node1, size_t node2, string &A)
{
    size_t inlet = 0, outlet = 0;
    if (this->dir_ == 0)
        inlet = 0,
        outlet = 1;
    else if (this->dir_ == 1)
        inlet = 2,
        outlet = 3;
    else if (this->dir_ == 2)
        inlet = 5,
        outlet = 4;

    if (_Pnt_attri[node1 - 1][inlet] == true &&
        _Pnt_attri[node2 - 1][inlet] == true)
    {
        A = "in";
        return true;
    }

    if (_Pnt_attri[node1 - 1][outlet] == true &&
        _Pnt_attri[node2 - 1][outlet] == true)
    {
        A = "out";
        return true;
    }

    return false;
};

inline bool Mesh_DFN_linear::If_interior_edge(size_t node1, size_t node2)
{
    bool A = true, B = true; // both are interior points

    for (size_t i = 0; i < 7; ++i)
        if (_Pnt_attri[node1 - 1][i] == true)
        {
            A = false;
            break;
        }

    for (size_t i = 0; i < 7; ++i)
        if (_Pnt_attri[node2 - 1][i] == true)
        {
            B = false;
            break;
        }

    if (A == true || B == true)
        return true;
    else
        return false;
};

inline void Mesh_DFN_linear::Identify_point_attribute(DFN::Domain dom)
{

    _Pnt_attri.resize(this->coordinate_3D.rows());
    for (size_t i = 0; i < (size_t)this->coordinate_3D.rows(); ++i)
    {
        _Pnt_attri[i] << false, false, false, false, false, false, false;

        Vector3d thisPnT;
        thisPnT << this->coordinate_3D(i, 0),
            this->coordinate_3D(i, 1),
            this->coordinate_3D(i, 2);

        // model bounds
        double top_zmax = dom.Model_domain(0);
        double bottom_zmin = dom.Model_domain(1);
        double front_ymin = dom.Model_domain(2);
        double back_ymax = dom.Model_domain(3);
        double left_xmin = dom.Model_domain(4);
        double right_xmax = dom.Model_domain(5);

        if (abs(thisPnT(2) - top_zmax) < 1e-4)
            _Pnt_attri[i][5] = true;

        if (abs(thisPnT(2) - bottom_zmin) < 1e-4)
            _Pnt_attri[i][4] = true;

        if (abs(thisPnT(1) - front_ymin) < 1e-4)
            _Pnt_attri[i][3] = true;

        if (abs(thisPnT(1) - back_ymax) < 1e-4)
            _Pnt_attri[i][2] = true;

        if (abs(thisPnT(0) - left_xmin) < 1e-4)
            _Pnt_attri[i][1] = true;

        if (abs(thisPnT(0) - right_xmax) < 1e-4)
            _Pnt_attri[i][0] = true;

        DFN::Point_3D ThisPnt{thisPnT};
        for (size_t k = 0; k < this->Frac_Tag.size(); ++k)
        {
            size_t FracID_ = Frac_Tag[k];
            if (ThisPnt.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true)
            {
                _Pnt_attri[i][6] = true;
                break;
            };
        }
    };
};

inline void Mesh_DFN_linear::Identify_ele_fracNO(DFN::Domain dom)
{
    this->element_2D.resize(this->Frac_Tag.size());

    vector<size_t> Num_ele_in_frac(this->Frac_Tag.size(), 0);

    vector<pair<size_t, size_t>> Tag__ele__(this->element_3D.size());

    //---------
    size_t Nproc = this->Nproc_global;
#pragma omp parallel for schedule(dynamic) num_threads(Nproc)
    for (int i = 0; i < this->element_3D.rows(); ++i)
    {
        size_t node1 = this->element_3D(i, 0),
               node2 = this->element_3D(i, 1),
               node3 = this->element_3D(i, 2);

        Vector3d A, B, C;
        A << (double)coordinate_3D(node1 - 1, 0),
            (double)coordinate_3D(node1 - 1, 1),
            (double)coordinate_3D(node1 - 1, 2);

        B << (double)coordinate_3D(node2 - 1, 0),
            (double)coordinate_3D(node2 - 1, 1),
            (double)coordinate_3D(node2 - 1, 2);

        C << (double)coordinate_3D(node3 - 1, 0),
            (double)coordinate_3D(node3 - 1, 1),
            (double)coordinate_3D(node3 - 1, 2);

        DFN::Point_3D ElePnt_1{A};
        DFN::Point_3D ElePnt_2{B};
        DFN::Point_3D ElePnt_3{C};

        for (size_t k = 0; k < this->Frac_Tag.size(); ++k)
        {
            size_t FracID_ = this->Frac_Tag[k];
            if ((ElePnt_1.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true || ElePnt_1.If_lies_within_a_polygon_3D(dom.Fractures[FracID_].Verts_trim) == true) &&
                (ElePnt_2.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true || ElePnt_2.If_lies_within_a_polygon_3D(dom.Fractures[FracID_].Verts_trim) == true) &&
                (ElePnt_3.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true || ElePnt_3.If_lies_within_a_polygon_3D(dom.Fractures[FracID_].Verts_trim) == true))
            {
                {
#pragma omp critical
                    {
                        Num_ele_in_frac[k]++;
                        Tag__ele__[i].second = Num_ele_in_frac[k] - 1;
                    }
                    Tag__ele__[i].first = k;
                }
                break;
            }

            if (k == this->Frac_Tag.size() - 1)
                throw Error_throw_ignore("Cannot find which frac does this element lie on. In class Mesh_DFN_overall");
        }
    }

#pragma omp parallel for schedule(dynamic) num_threads(Nproc)
    for (size_t i = 0; i < this->Frac_Tag.size(); ++i)
        this->element_2D[i].resize(Num_ele_in_frac[i], 3);

#pragma omp parallel for schedule(dynamic) num_threads(Nproc)
    for (int i = 0; i < this->element_3D.rows(); ++i)
        this->element_2D[Tag__ele__[i].first].row(Tag__ele__[i].second) = this->element_3D.row(i);
};

inline SparseMatrix<double> Mesh_DFN_linear::Node2element(size_t i)
{
    size_t Dim = this->coordinate_3D.rows();
    SparseMatrix<double> Node2element__k(Dim, Dim);
    Node2element__k.reserve(VectorXi::Constant(Dim, 5));

    for (int j = 0; j < this->element_2D[i].rows(); ++j)
    {
        for (size_t k = 0; k < 3; ++k)
        {
            size_t node1 = this->element_2D[i](j, k);
            size_t node2 = this->element_2D[i](j, (k + 1) % 3);

            /*
            if (Node2element__k.coeffRef(node1 - 1, node2 - 1) != 0)
            {
                cout << node1 << ", " << node2 << endl;

                cout << this->element_2D[i] << endl;
                this->Matlab_plot("YTC.mat", "YTC.m", dom);
                exit(0);
            }*/

            Node2element__k.coeffRef(node1 - 1, node2 - 1) = j + 1;
        }
    }

    Node2element__k.makeCompressed();

    return Node2element__k;
};

void Mesh_DFN_linear::Matlab_plot(string FileKey_mat,
                                  string FileKey_m,
                                  DFN::Domain dom)
{
    const char *filename = FileKey_mat.c_str();

    DFN::MATLAB_DATA_API M1_;
    M1_.Write_mat(filename, "w", 1, 1, 1, {0}, "nothing_");

    for (size_t i = 0; i < this->Frac_Tag.size(); ++i)
    {
        //cout << i << endl;
        MatrixXs ele_2D_frac = this->element_2D[i];

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
        /*
        vector<double> pData2(this->coordinate_2D[i].rows() * 2);

        for (size_t j = 0; j < (size_t)this->coordinate_2D[i].rows() * 2; ++j)
        {
            size_t k, l;
            k = ceil(j / this->coordinate_2D[i].rows()); // column
            l = j % this->coordinate_2D[i].rows();       // row

            pData2[j] = this->coordinate_2D[i].coeffRef(l, k);
        }
        M1_.Write_mat(filename, "u", this->coordinate_2D[i].rows() * 2,
                      this->coordinate_2D[i].rows(), 2, pData2,
                      "coordinate_2D_Frac_" + to_string(i + 1));*/

        //----------------------------------
    }

    vector<double> pData4(this->coordinate_3D.rows() * 3);
    vector<double> pData5(this->element_3D.rows() * 3);

    for (size_t j = 0; j < (size_t)this->coordinate_3D.rows() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / this->coordinate_3D.rows()); // column
        l = j % this->coordinate_3D.rows();       // row

        pData4[j] = this->coordinate_3D(l, k);
        //cout << pData4[j] << endl;
    }

    for (size_t j = 0; j < (size_t)this->element_3D.rows() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / this->element_3D.rows()); // column
        l = j % this->element_3D.rows();       // row

        pData5[j] = this->element_3D(l, k);
    }

    string FracJXY3D_s = "coordinate_3D";
    string FracJM_s = "element_3D";

    M1_.Write_mat(filename, "u", this->coordinate_3D.rows() * 3,
                  this->coordinate_3D.rows(), 3, pData4, FracJXY3D_s);
    M1_.Write_mat(filename, "u", this->element_3D.rows() * 3,
                  this->element_3D.rows(), 3, pData5, FracJM_s);

    //---------------------
    std::map<pair<size_t, size_t>, size_t> Inlet_edges, Outlet_edges, Neumann_edges, Interior_edgeNO;
    for (size_t i = 0; i < element_2D.size(); ++i)
        for (int j = 0; j < element_2D[i].rows(); ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                size_t node1 = element_2D[i](j, k),
                       node2 = element_2D[i](j, (k + 1) % 3);
                std::pair<size_t, size_t> At_ = std::make_pair(node1 < node2 ? node1 : node2,
                                                               node1 > node2 ? node1 : node2);

                if (this->Edge_attri_and_NO[i](j, k).first == "in")
                    Inlet_edges.insert(std::make_pair(At_, 0));
                else if (this->Edge_attri_and_NO[i](j, k).first == "out")
                    Outlet_edges.insert(std::make_pair(At_, 0));
                else if (this->Edge_attri_and_NO[i](j, k).first == "neumann")
                    Neumann_edges.insert(std::make_pair(At_, 0));
                else
                    Interior_edgeNO.insert(std::make_pair(At_, 0));
            }
        }

    vector<double> pData6(Inlet_edges.size() * 2);
    size_t yw = 0;
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Inlet_edges.begin();
         its != Inlet_edges.end(); its++)
    {
        size_t node1 = its->first.first;
        pData6[yw] = node1;
        yw++;
    }
    yw = Inlet_edges.size();
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Inlet_edges.begin();
         its != Inlet_edges.end(); its++)
    {
        size_t node1 = its->first.second;
        pData6[yw] = node1;
        yw++;
    }

    M1_.Write_mat(filename, "u", Inlet_edges.size() * 2,
                  Inlet_edges.size(), 2, pData6, "Inlet_edges");

    //-------------------------
    vector<double> pData7(Outlet_edges.size() * 2);
    yw = 0;
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Outlet_edges.begin();
         its != Outlet_edges.end(); its++)
    {
        size_t node1 = its->first.first;
        pData7[yw] = node1;
        yw++;
    }
    yw = Outlet_edges.size();
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Outlet_edges.begin();
         its != Outlet_edges.end(); its++)
    {
        size_t node1 = its->first.second;
        pData7[yw] = node1;
        yw++;
    }

    M1_.Write_mat(filename, "u", Outlet_edges.size() * 2,
                  Outlet_edges.size(), 2, pData7, "Outlet_edges");

    //-------------------------
    vector<double> pData8(Neumann_edges.size() * 2);
    yw = 0;
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Neumann_edges.begin();
         its != Neumann_edges.end(); its++)
    {
        size_t node1 = its->first.first;
        pData8[yw] = node1;
        yw++;
    }
    yw = Neumann_edges.size();
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Neumann_edges.begin();
         its != Neumann_edges.end(); its++)
    {
        size_t node1 = its->first.second;
        pData8[yw] = node1;
        yw++;
    }

    M1_.Write_mat(filename, "u", Neumann_edges.size() * 2,
                  Neumann_edges.size(), 2, pData8, "Neumann_edges");

    vector<double> pData9(Interior_edgeNO.size() * 2);
    yw = 0;
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Interior_edgeNO.begin();
         its != Interior_edgeNO.end(); its++)
    {
        size_t node1 = its->first.first;
        pData9[yw] = node1;
        yw++;
    }
    yw = Interior_edgeNO.size();
    for (std::map<pair<size_t, size_t>, size_t>::iterator its = Interior_edgeNO.begin();
         its != Interior_edgeNO.end(); its++)
    {
        size_t node1 = its->first.second;
        pData9[yw] = node1;
        yw++;
    }

    M1_.Write_mat(filename, "u", Interior_edgeNO.size() * 2,
                  Interior_edgeNO.size(), 2, pData9, "Interior_edgeNO");

    // m file
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;";
    oss << "load('" << FileKey_mat << "');\n";

    oss << "figure(1); patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', zeros(size(coordinate_3D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);\n";
    oss << "view(3); hold on;\n";
    oss << "for i = 1:size(coordinate_3D, 1)\n";
    oss << "\t;\n";
    oss << "end\n";
    oss << "% text(coordinate_3D(:, 1), coordinate_3D(:, 2), coordinate_3D(:, 3), num2str([1:1:size(coordinate_3D, 1)]')); hold on\n";
    oss << "title('mesh result')\n\n";

    oss << "figure(2)\n";
    for (size_t i = 0; i < this->Frac_Tag.size(); i++)
    {
        size_t nf = this->Frac_Tag[i];
        size_t n_verts = dom.Fractures[nf].Verts_trim.size();
        oss << "frac" << nf + 1 << " = fill3([";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << dom.Fractures[nf].Verts_trim[nv_1](0) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << dom.Fractures[nf].Verts_trim[nv_1](1) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << dom.Fractures[nf].Verts_trim[nv_1](2) << " ";
        }
        oss << "],[rand rand rand]);\ngrid on;\nhold on;\n";
    }

    oss << "\nA = []; B = []; C = [];\n";
    oss << "for i = 1:size(Inlet_edges, 1)\n";
    oss << "\tnode1 = Inlet_edges(i, 1); node2 = Inlet_edges(i, 2);\n";
    oss << "\tA(:, i) = [coordinate_3D(node1, 1), coordinate_3D(node2, 1)]';\n";
    oss << "\tB(:, i) = [coordinate_3D(node1, 2), coordinate_3D(node2, 2)]';\n";
    oss << "\tC(:, i) = [coordinate_3D(node1, 3), coordinate_3D(node2, 3)]';\n";
    oss << "end\n";
    oss << "plot3(A, B, C, 'r-*', 'linewidth', 1.5); hold on\n";

    oss << "\nA = []; B = []; C = [];\n";
    oss << "for i = 1:size(Outlet_edges, 1)\n";
    oss << "\tnode1 = Outlet_edges(i, 1); node2 = Outlet_edges(i, 2);\n";
    oss << "\tA(:, i) = [coordinate_3D(node1, 1), coordinate_3D(node2, 1)]';\n";
    oss << "\tB(:, i) = [coordinate_3D(node1, 2), coordinate_3D(node2, 2)]';\n";
    oss << "\tC(:, i) = [coordinate_3D(node1, 3), coordinate_3D(node2, 3)]';\n";
    oss << "end\n";
    oss << "plot3(A, B, C, 'b-o', 'linewidth', 1.5); hold on\n";

    oss << "\nA = []; B = []; C = [];\n";
    oss << "for i = 1:size(Neumann_edges, 1)\n";
    oss << "\tnode1 = Neumann_edges(i, 1); node2 = Neumann_edges(i, 2);\n";
    oss << "\tA(:, i) = [coordinate_3D(node1, 1), coordinate_3D(node2, 1)]';\n";
    oss << "\tB(:, i) = [coordinate_3D(node1, 2), coordinate_3D(node2, 2)]';\n";
    oss << "\tC(:, i) = [coordinate_3D(node1, 3), coordinate_3D(node2, 3)]';\n";
    oss << "end\n";
    oss << "% plot3(A, B, C, 'k-', 'linewidth', 1.5); hold on\n";

    oss << "\nA = []; B = []; C = [];\n";
    oss << "for i = 1:size(Interior_edgeNO, 1)\n";
    oss << "\tnode1 = Interior_edgeNO(i, 1); node2 = Interior_edgeNO(i, 2);\n";
    oss << "\tA(:, i) = [coordinate_3D(node1, 1), coordinate_3D(node2, 1)]';\n";
    oss << "\tB(:, i) = [coordinate_3D(node1, 2), coordinate_3D(node2, 2)]';\n";
    oss << "\tC(:, i) = [coordinate_3D(node1, 3), coordinate_3D(node2, 3)]';\n";
    oss << "end\n";
    oss << "% plot3(A, B, C, 'g-', 'linewidth', 1.5); hold on\n";

    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << dom.Surfaces[i].Verts[j](0) << " " << dom.Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](1) << " " << dom.Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](2) << " " << dom.Surfaces[i].Verts[nj](2) << "],";
            oss << "'k-','Linewidth',3);\ngrid on; hold on;\n";
        }
    }
    double xmin_1 = dom.Model_domain(4), xmax_1 = dom.Model_domain(5);
    double ymin_1 = dom.Model_domain(2), ymax_1 = dom.Model_domain(3);
    double zmin_1 = dom.Model_domain(1), zmax_1 = dom.Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

    oss.close();
};

inline size_t Mesh_DFN_linear::Element_ID(size_t i, size_t j)
{
    size_t ID = 0;

    if (i == 0)
        return ID + j;

    for (size_t y = 0; y < i; ++y)
        ID += this->element_2D[y].rows();

    return ID + j;
};
} // namespace DFN
