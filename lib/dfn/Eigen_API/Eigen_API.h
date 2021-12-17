#pragma once
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <algorithm>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

typedef Matrix<size_t, Dynamic, Dynamic> MatrixXs;
typedef SparseMatrix<float> Sp_f_mat;
typedef SparseMatrix<size_t> Sp_s_mat;
typedef Matrix<bool, 7, 1> Vector7b;

MatrixXs removeMatrixXsRow(const MatrixXs original_matrix,
                           vector<int> row_to_remove);

//--------------

MatrixXs removeMatrixXsRow(const MatrixXs original_matrix,
                           vector<int> row_to_remove)
{
    // New matrix has one fewer rows
    MatrixXs new_matrix(original_matrix.rows() - row_to_remove.size(), original_matrix.cols());
    // Track rows in new matrix. Skip one at row_to_remove.
    int row_to_fill = 0;
    
    for (int orig_matrix_row = 0;
         orig_matrix_row < original_matrix.rows();
         ++orig_matrix_row)
    {
        vector<int>::iterator its =
            find(row_to_remove.begin(), row_to_remove.end(), orig_matrix_row);

        if (its == row_to_remove.end())
        {
            new_matrix.row(row_to_fill) = original_matrix.row(orig_matrix_row);
            ++row_to_fill;
        }
    }
    return new_matrix;
}