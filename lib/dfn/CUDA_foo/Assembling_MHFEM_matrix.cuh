#ifndef ASSEMBLINE_MHFEM_MATRIX_CUH
#define ASSEMBLINE_MHFEM_MATRIX_CUH

#include "../Mesh_H/Mesh_DFN_linear.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

extern "C"
{
    Eigen::SparseMatrix<double> Assembling_MHFEM_matrix(DFN::Mesh_DFN_linear mesh,
                                                        DFN::Domain dom);
}

#endif