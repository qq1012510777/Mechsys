#pragma once
#include "../Mesh_DFN_linear.h"

namespace DFN
{
class GLOB_edge_numbering
{
public:
    vector<MatrixXs> GLOB_edges_NO_of_ele_frac;

    Sp_s_mat nodes2GLOB_edge;
    Sp_s_mat nodes2GLOB_interior_edge;
    size_t NUM_GLOB_edge = 0;
    size_t NUM_interior_edges = 0;

public:
    GLOB_edge_numbering(DFN::Mesh_DFN_linear mesh);
};

inline GLOB_edge_numbering::GLOB_edge_numbering(DFN::Mesh_DFN_linear mesh)
{

    GLOB_edges_NO_of_ele_frac.resize(mesh.element_2D.size());

    SparseMatrix<size_t> m1(mesh.coordinate_3D.rows(),
                            mesh.coordinate_3D.rows());
    m1.reserve(VectorXi::Constant((int)mesh.coordinate_3D.rows(), 4));

    size_t Global_edgeNO = 1;

    for (size_t i = 0; i < (size_t)mesh.element_3D.rows(); ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            size_t node1 = mesh.element_3D(i, j),
                   node2 = mesh.element_3D(i, (j + 1) % 3);

            if (m1.coeffRef(node1 - 1, node2 - 1) == 0 &&
                m1.coeffRef(node2 - 1, node1 - 1) == 0)
            {
                m1.coeffRef(node1 - 1, node2 - 1) = Global_edgeNO;
                m1.coeffRef(node2 - 1, node1 - 1) = Global_edgeNO;
                Global_edgeNO++;
            }
        }
    }
    m1.makeCompressed();
    nodes2GLOB_edge = m1;

    NUM_GLOB_edge = Global_edgeNO - 1;

    for (size_t i = 0; i < GLOB_edges_NO_of_ele_frac.size(); ++i)
    {
        GLOB_edges_NO_of_ele_frac[i].resize(mesh.element_2D[i].rows(),
                                            mesh.element_2D[i].cols());

        for (size_t j = 0; j < (size_t)mesh.element_2D[i].rows(); ++j)
            for (size_t k = 0; k < 3; ++k)
            {
                size_t node1 = mesh.element_2D[i](j, k),
                       node2 = mesh.element_2D[i](j, (k + 1) % 3);

                size_t GLOB_edgeNO = this->nodes2GLOB_edge.coeffRef(node1 - 1,
                                                                    node2 - 1);

                GLOB_edges_NO_of_ele_frac[i](j, k) = GLOB_edgeNO;
            }
    };

    SparseMatrix<size_t> m2(mesh.coordinate_3D.rows(),
                            mesh.coordinate_3D.rows());
    m2.reserve(VectorXi::Constant((int)mesh.coordinate_3D.rows(), 4));

    for (size_t i = 0; i < (size_t)mesh.Interior_edges.size(); ++i)
    {
        for (size_t j = 0; j < (size_t)mesh.Interior_edges[i].rows(); ++j)
        {
            size_t node1 = mesh.Interior_edges[i](j, 0),
                   node2 = mesh.Interior_edges[i](j, 1);

            if (m2.coeffRef(node1 - 1, node2 - 1) == 0)
            {
                NUM_interior_edges++;
                m2.coeffRef(node1 - 1, node2 - 1) = NUM_interior_edges;
                m2.coeffRef(node2 - 1, node1 - 1) = NUM_interior_edges;
            }
        }
    }

    m2.makeCompressed();
    nodes2GLOB_interior_edge = m2;
};
}; // namespace DFN