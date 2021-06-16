#pragma once
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
    double *X_K;
public:
    Using_UMFPACK(const double *K_a,
                  const size_t Dim /*dimension of the square matrix*/,
                  const double *B);
    ~Using_UMFPACK();
};

inline Using_UMFPACK::Using_UMFPACK(const double *K_a, const size_t Dim, const double *B)
{
    std::vector<double> Ax_xx;

    std::vector<int> Ai_xx;

    std::vector<int> Ap_xx;

    for (size_t is = 0; is < Dim; ++is)
    {
        if (is == 0)
            Ap_xx.push_back(0);
        else
            Ap_xx.push_back(Ai_xx.size());

        for (size_t js = 0; js < Dim; ++js)
        {
            size_t k = js * Dim + is; // column major
            if (abs(K_a[k]) > 1e-6)
            {
                Ax_xx.push_back(K_a[k]);
                Ai_xx.push_back(js);
            }
        }
    }

    Ap_xx.push_back(Ai_xx.size());

    int *Ai;
    int *Ap;
    double *Ax;

    Ai = (int *)calloc(Ai_xx.size(), sizeof(int));
    Ap = (int *)calloc(Ap_xx.size(), sizeof(int));
    Ax = (double *)calloc(Ax_xx.size(), sizeof(double));

    for (size_t is = 0; is < Ai_xx.size(); ++is)
    {
        Ai[is] = Ai_xx[is];
        Ax[is] = Ax_xx[is];
    }

    for (size_t is = 0; is < Ap_xx.size(); ++is)
    {
        Ap[is] = Ap_xx[is];
    }

    int n = Dim;
    double *null = (double *)NULL;
    void *Numeric;
    int status;
    void *Symbolic;
    double x[Dim];

    status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, null, null);
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic(&Symbolic);
    status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, B, Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);

    X_K = (double *)calloc(Dim, sizeof(double));

    for (size_t is = 0; is < Dim; is++)
        X_K[is] = x[is];   

    free(Ai);
    free(Ap);
    free(Ax);

    status++;
};

inline Using_UMFPACK::~Using_UMFPACK()
{
    //cout << "free X_K\n";
    free(X_K);
};

} // namespace DFN