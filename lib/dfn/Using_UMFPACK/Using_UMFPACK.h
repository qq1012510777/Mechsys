#pragma once
#include "../Error_throw/Error_throw.h"
#include "amd.h"
#include "umfpack.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

namespace DFN
{
class Using_UMFPACK
{
public:
    int *Ai;
    int *Ap;
    double *Ax;

public:
    Using_UMFPACK();

    void Prepare(const float *K_a, const size_t Dim);

    void Solve(/*dimension of the square matrix*/
               const size_t Dim,
               double *B);
};

inline Using_UMFPACK::Using_UMFPACK(){

};

inline void Using_UMFPACK::Prepare(const float *K_a, const size_t Dim)
{
    size_t Size_AP = 0, Size_Ax = 0, Size_Ai = 0;

    Size_AP = Dim + 1;

    for (size_t is = 0; is < Dim; ++is)
    {
        for (size_t js = 0; js < Dim; ++js)
        {
            size_t k = js * Dim + is;

            if (abs(K_a[k]) > 1e-6)
            {
                Size_Ax++;
                Size_Ai++;
            }
        }
    }

    //cout << "umfpack mat init start_\n";
    Ai = new int[Size_Ai]();

    if (Ai == NULL)
    {
        cout << "\t\tNUll Ai\n";
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'Ai'!\n");
    }

    Ap = new int[Size_AP]();

    if (Ap == NULL)
    {
        cout << "\t\tNUll Ap\n";
        delete[] Ai;
        Ai = NULL;
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'Ap'!\n");
    }

    Ax = new double[Size_Ax]();

    if (Ax == NULL)
    {
        delete[] Ai;
        Ai = NULL;
        delete[] Ap;
        Ap = NULL;
        cout << "\t\tNUll Ax\n";
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'Ax'!\n");
    }

    //cout << "umfpack mat init finish_\n";

    size_t Id_non_zero = 0;
    for (size_t is = 0; is < Dim; ++is)
    {
        if (is == 0)
            Ap[is] = 0;
        else
            Ap[is] = Id_non_zero;

        for (size_t js = 0; js < Dim; ++js)
        {
            size_t k = js * Dim + is; // column major
            if (abs(K_a[k]) > 1e-6)
            {
                Ax[Id_non_zero] = (double)K_a[k];
                Ai[Id_non_zero] = (int)js;
                Id_non_zero++;
            }
        }
    }
    Ap[Size_AP - 1] = Size_Ax;
};

inline void Using_UMFPACK::Solve(const size_t Dim, double *B)
{
    int n = Dim;
    double *null = (double *)NULL;
    void *Numeric;
    int status;
    void *Symbolic;

    double *x = new double[Dim]();

    if (x == NULL)
    {
        delete[] Ai;
        Ai = NULL;
        delete[] Ap;
        Ap = NULL;
        delete[] Ax;
        Ax = NULL;
        cout << "\t\tNUll x\n";
        throw Error_throw_ignore("Error! Cannot alloc to matrix 'x'!\n");
    }

    status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, null, null);
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic(&Symbolic);
    status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, B, Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);

    delete[] Ai;
    Ai = NULL;
    delete[] Ap;
    Ap = NULL;
    delete[] Ax;
    Ax = NULL;

    for (size_t is = 0; is < Dim; is++)
        B[is] = x[is];

    delete[] x;
    x = NULL;

    status++;
};

} // namespace DFN