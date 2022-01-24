#pragma once
#include "../Error_throw/Error_throw.h"
#include "amd.h"
#include "umfpack.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{
class Using_UMFPACK_Eigen
{
public:
    Using_UMFPACK_Eigen(SparseMatrix<double> D,
                        SparseMatrix<double> r,
                        MatrixXd &p);
};

inline Using_UMFPACK_Eigen::Using_UMFPACK_Eigen(SparseMatrix<double> D,
                                                SparseMatrix<double> r,
                                                MatrixXd &p)
{
    int n = D.rows();
    double *null = (double *)NULL;
    void *Numeric;
    int status;
    void *Symbolic;
    //
    //
    size_t NUM_nonZeros = D.nonZeros();

    int *Ai = new int[NUM_nonZeros]();
    
    if (Ai == NULL)
    {
        //delete[] Ai;
        //Ai = NULL;
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'x'!\n");
    }

    int *Ap = new int[n + 1]();
    if (Ap == NULL)
    {
        delete[] Ai;
        Ai = NULL;
        //delete[] Ap;
        //Ap = NULL;
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'Ap'!\n");
    }

    double *Ax = new double[NUM_nonZeros]();
    if (Ax == NULL)
    {
        delete[] Ai;
        Ai = NULL;
        delete[] Ap;
        Ap = NULL;
        //delete[] Ax;
        //Ax = NULL;
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'Ax'!\n");
    }

    double *b = new double[n]();
    if (b == NULL)
    {
        delete[] Ai;
        Ai = NULL;
        delete[] Ap;
        Ap = NULL;
        delete[] Ax;
        Ax = NULL;
        //delete[] b;
        //b = NULL;
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'b'!\n");
    }

    double *x = new double[n]();
    if (x == NULL)
    {
        delete[] Ai;
        Ai = NULL;
        delete[] Ap;
        Ap = NULL;
        delete[] Ax;
        Ax = NULL;
        delete[] b;
        b = NULL;
        //delete[] x;
        //x = NULL;
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'x'!\n");
    }

    VectorXd RT = VectorXd(r);
    std::copy(RT.data(), RT.data() + RT.size(), b);
    //cout << 1 << endl;
    size_t tmp_Id = 0;
    for (int k = 0; k < D.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(D, k); it; ++it)
        {
            Ai[tmp_Id] = it.row();
            Ax[tmp_Id] = it.value();
            tmp_Id++;
        }
        Ap[k + 1] = tmp_Id;
    }

    status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, null, null);
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic(&Symbolic);
    status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);
    //cout << 1.1 << endl;
    p.resize(n, 1);
    p << Map<VectorXd>(x, n);
    //cout << 2 << endl;
    delete[] Ai;
    Ai = NULL;
    delete[] Ap;
    Ap = NULL;
    delete[] Ax;
    Ax = NULL;
    delete[] b;
    b = NULL;
    delete[] x;
    x = NULL;

    status++;
};

}; // namespace DFN