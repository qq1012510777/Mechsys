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
    Using_UMFPACK(const double *K_a,
                  const size_t Dim /*dimension of the square matrix*/,
                  double *B);
};

inline Using_UMFPACK::Using_UMFPACK(const double *K_a, const size_t Dim, double *B)
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

    int *Ai = Ai_xx.data();
    int *Ap = Ap_xx.data();
    double *Ax = Ax_xx.data();

    Ai_xx.clear();
    Ap_xx.clear();
    Ax_xx.clear();

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

    for (size_t is = 0; is < Dim; is++)
        B[is] = x[is];

    status++;
};

} // namespace DFN