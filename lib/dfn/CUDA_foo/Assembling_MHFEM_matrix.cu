#include "Assembling_MHFEM_matrix.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void Assemble_on_GPU();

Eigen::SparseMatrix<double> Assembling_MHFEM_matrix(DFN::Mesh_DFN_linear mesh,
                                                    DFN::Domain dom)
{
    //
    Eigen::SparseMatrix<double> K;

    return K;
}