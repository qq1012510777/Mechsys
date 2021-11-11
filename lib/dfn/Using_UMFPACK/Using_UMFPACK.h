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
    std::vector<double> Ax_xx;

    std::vector<int> Ai_xx;

    std::vector<int> Ap_xx;

public:
    Using_UMFPACK();

    void Prepare(const double *K_a,
                 const size_t Dim);

    void Solve(const size_t Dim /*dimension of the square matrix*/,
               double *B);
};

inline Using_UMFPACK::Using_UMFPACK()
{
    Ax_xx.resize(0);
    Ai_xx.resize(0);
    Ap_xx.resize(0);
};

inline void Using_UMFPACK::Prepare(const double *K_a,
                                   const size_t Dim)
{
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
}

inline void Using_UMFPACK::Solve(const size_t Dim, double *B)
{

    //cout << "umfpack mat init start_\n";
    int *Ai = (int *)malloc(Ai_xx.size() * sizeof(int));
    memset(Ai, 0, Ai_xx.size() * sizeof(int));
    int *Ap = (int *)malloc(Ap_xx.size() * sizeof(int));
    memset(Ap, 0, Ap_xx.size() * sizeof(int));
    double *Ax = (double *)malloc(Ax_xx.size() * sizeof(double));
    memset(Ax, 0, Ax_xx.size() * sizeof(double));
    //cout << "umfpack mat init finish_\n";

    std::copy(Ai_xx.begin(), Ai_xx.end(), Ai);
    std::copy(Ap_xx.begin(), Ap_xx.end(), Ap);
    std::copy(Ax_xx.begin(), Ax_xx.end(), Ax);

    this->Ai_xx.clear();
    this->Ap_xx.clear();
    this->Ax_xx.clear();

    int n = Dim;
    double *null = (double *)NULL;
    void *Numeric;
    int status;
    void *Symbolic;

    double *x = (double *)malloc(Dim * sizeof(double));
    memset(x, 0, Dim * sizeof(double));

    status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, null, null);
    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);
    umfpack_di_free_symbolic(&Symbolic);
    status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, B, Numeric, null, null);
    umfpack_di_free_numeric(&Numeric);

    free(Ai);
    free(Ap);
    free(Ax);

    for (size_t is = 0; is < Dim; is++)
        B[is] = x[is];

    free(x);

    status++;
};

} // namespace DFN