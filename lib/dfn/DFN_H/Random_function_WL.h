#pragma once
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include "../Error_throw/Error_throw.h"

namespace DFN
{
using namespace std;
/*------------
part of this code is copied from GeneralBlock
http://www.rockfractures.com/
Dr. Qingchun Yu (China University of Geosciences (Beijing))
*/
class Random_function
{
private:
    gsl_rng *dmodul;

public:
    Random_function(gsl_rng *seed_factor);                             ///< change the "seed_factor" to generate different data series
                                                                       ///< " 0.0 < seed_factor < 1.0"
    double rn();                                                       ///< 0-1 uniform random data
    double unifrm(double a, double b);                                 ///< a-b uniform random data
    double gauss(double mean, double sigma, double min, double max);   ///< normal distributed random data
    double lognor(double mean, double sigma, double min, double max);  ///< log_normal distributed data
    double erlang(double mean, double min, double max);                ///< exponential distributed data
    void fisher(double mean_y, double mean_x, double fisher_k, ///< Fisher, 1st parameter: dip direction, 2nd: dip angle, 3rd: fisher constant,
                        double min_dip_dir, double max_dip_dir,
                        double min_dip_ang, double max_dip_ang,
                        double &dd_alpha1, double &da_beta1);

    void lubksb(int n, int np, int *indx, double *b);           ///< the later 2 are used only when transform (0,0) to (meanx,meany)
    void ludcmp(int n, int np, int *indx);
    double fisher_a[10][10]; ///< a[][] was used only in the fisher-generator

    double powerlaw(double x0, double x1, double alpha_g); ///< powerlaw
};

inline Random_function::Random_function(gsl_rng *seed_factor)
{
    dmodul = seed_factor;
    //cout << "random seed: " << dmodul << endl;
    return;
}
//---------------------------------------------
inline double Random_function::rn()
{

    double rrn = gsl_rng_uniform(dmodul);
    return rrn;
}
//-------------------- ---------------------------
inline double Random_function::unifrm(double a, double b)
{
    double random, rrn;
    rrn = rn();
    //std::cout << rrn<<"\n";
    random = a + (b - a) * rrn;
    return random;
}
//------------------------------------------
inline double Random_function::gauss(double mean, double sigma, double min, double max)
{
    double v, r1, r2, random;
g100:;
    r1 = rn(); //
    r2 = rn();
    v = sqrt((-2. * log(r1))) * cos(6.2832 * r2);
    random = v * sigma + mean;
    if (random < min || random > max)
        goto g100;
    return random;
}
//----------------------------------------------
inline double Random_function::lognor(double mean, double sigma, double min, double max) //sigma is variance
{

    double random;

    double mean_1 = log(mean * mean / (pow(sigma + mean * mean, 0.5)));
    double sigma_1 = pow(log(1 + ((double)sigma) / (mean * mean)), 0.5); //sigma_1 is input std. deviation
l100:;
    random = gsl_ran_lognormal(dmodul, mean_1, sigma_1);
    if (random < min || random > max)
        goto l100;
    return random;
}
//-------------------------------------------------
inline double Random_function::erlang(double mean, double min, double max)
{
    double r, random, alpha;
    alpha = 1. / mean;
e100:;
    r = 1.0;
    r = r * rn();
    if (r == 0.)
        goto e100;
    random = -1.0 / alpha * log(r);
    if (random < min || random > max)
        goto e100;
    return random;
}
//----------------------------------------------
inline void Random_function::fisher(double mean_y, double mean_x, double fisher_k,
                                            double min_dip_dir, double max_dip_dir,
                                            double min_dip_ang, double max_dip_ang,
                                            double &dd_alpha1, double &da_beta1)
{
    int n;
    double dip_dir1;
    double dip_ang1;
    double rrn, PI, fisher_x, fisher_y;
    int indx[10];
    double b[10];
    double l1, m1, n1, lo, mo, no, phi, seta, dip_dir, dip_ang;
    //printf("------\n%f, %f, %f, %f, %f, %f, %f\n", mean_y, mean_x, fisher_k, min_dip_dir, max_dip_dir,min_dip_ang, max_dip_ang);
    PI = 3.14159;
    mean_y = mean_y / 180.0 * PI;
    mean_x = mean_x / 180.0 * PI;
    mean_y = 2.5 * PI - mean_y;
    if (mean_y > 2.0 * PI)
        mean_y = mean_y - 2.0 * PI; //from dip direction to fisher-seta

    n = 0;
    while (n <= 0)
    {
        rrn = rn(); // A 0.0--1.0 random number
        //cout<<"rrn: "<<rrn<<endl;
        fisher_x = log(exp(fisher_k) - (exp(fisher_k) - exp(-fisher_k)) * rrn);
        fisher_x = fisher_x / fisher_k;
        fisher_x = acos(fisher_x);
        fisher_y = rn() * 2.0 * PI; // A random number 0.0--2PI

        l1 = sin(fisher_x) * cos(fisher_y);
        m1 = sin(fisher_x) * sin(fisher_y);
        n1 = cos(fisher_x);

        fisher_a[1][1] = cos(mean_x) * cos(mean_y);
        fisher_a[1][2] = cos(mean_x) * sin(mean_y);
        fisher_a[1][3] = -sin(mean_x);
        fisher_a[2][1] = -sin(mean_y);
        fisher_a[2][2] = cos(mean_y);
        fisher_a[2][3] = 0.0;
        fisher_a[3][1] = sin(mean_x) * cos(mean_y);
        fisher_a[3][2] = sin(mean_x) * sin(mean_y);
        fisher_a[3][3] = cos(mean_x);
        b[1] = l1;
        b[2] = m1;
        b[3] = n1;

        ludcmp(3, 3, indx);
        lubksb(3, 3, indx, b); //LU decomposition for linear equations
        lo = b[1];
        mo = b[2];
        no = b[3];
        phi = acos(no);
        seta = lo / sin(phi);
        seta = acos(seta);
        //cout<<"seta:::"<<seta<<endl;
        if (mo < 0.0)
            seta = 2 * PI - seta;
        if (phi > 0.5 * PI)
        {
            phi = PI - phi;
            seta = seta + PI;
        }
        dip_dir = 2.5 * PI - seta;
        if (dip_dir >= 2.0 * PI)
            dip_dir = dip_dir - 2.0 * PI;
        dip_ang = phi;
        dip_dir = dip_dir / PI * 180.0;
        dip_ang = dip_ang / PI * 180.0;

        if (dip_dir >= min_dip_dir && dip_dir <= max_dip_dir &&
            dip_ang >= min_dip_ang && dip_ang <= max_dip_ang)
        {
            n = n + 1;
            dip_dir1 = dip_dir;
            dip_ang1 = dip_ang;
        }
        //fprintf(seep,"%15lf %15lf\n",dip_dir,dip_ang);
    }
    
    dd_alpha1 = dip_dir1;
    da_beta1 = dip_ang1;
}


//--------------------------------------
inline void Random_function::lubksb(int n, int np, int *indx, double *b)
{
    int i, ii, j, ll;
    double sum;
    ii = 0;
    for (i = 1; i <= n; i++)
    { //--12
        ll = *(indx + i);
        sum = *(b + ll);
        *(b + ll) = *(b + i);
        if (ii != 0)
        {
            for (j = ii; j <= i - 1; j++)
                sum = sum - fisher_a[i][j] * (*(b + j));
        }
        else
        {
            if (sum != 0.)
                ii = i;
        }

        *(b + i) = sum;
    } //12    continue
    for (i = n; i >= 1; i--)
    {
        sum = *(b + i);
        for (j = i + 1; j <= n; j++)
            sum = sum - fisher_a[i][j] * (*(b + j));
        *(b + i) = sum / fisher_a[i][i];
    }

    return;
}
//-------------------------------------------------------
inline void Random_function::ludcmp(int n, int np, int *indx)
{
#define NMAX 500
#define TINY 1.0e-20

    int i, imax = 0, j, k, d;
    double aamax, dum, sum, vv[NMAX];

    d = 1;
    for (i = 1; i <= n; i++)
    { //do 12 i=1,n
        aamax = 0.0;
        for (j = 1; j <= n; j++)
        {
            if (fabs(fisher_a[i][j]) > aamax)
                aamax = fabs(fisher_a[i][j]);
        }
        vv[i] = 1. / aamax;
    } //12    continue
    for (j = 1; j <= n; j++)
    { //do 19 j=1,n
        for (i = 1; i <= j - 1; i++)
        {
            sum = fisher_a[i][j];
            for (k = 1; k <= i - 1; k++)
                sum = sum - fisher_a[i][k] * fisher_a[k][j];
            fisher_a[i][j] = sum;
        }
        aamax = 0.0;
        for (i = j; i <= n; i++)
        { //do 16 i=j,n
            sum = fisher_a[i][j];
            for (k = 1; k <= j - 1; k++)
                sum = sum - fisher_a[i][k] * fisher_a[k][j];
            fisher_a[i][j] = sum;
            dum = vv[i] * fabs(sum);
            if (dum >= aamax)
            {
                imax = i;
                aamax = dum;
            }
        } //16      continue
        if (j != imax)
        {
            for (k = 1; k <= n; k++)
            {
                dum = fisher_a[imax][k];
                fisher_a[imax][k] = fisher_a[j][k];
                fisher_a[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        *(indx + j) = imax;
        if (fisher_a[j][j] == 0.0)
            fisher_a[j][j] = TINY;
        if (j != n)
        {
            dum = 1.0 / fisher_a[j][j];
            for (i = j + 1; i <= n; i++)
                fisher_a[i][j] = fisher_a[i][j] * dum;
        }
    } //19    continue
    return;
}

inline double Random_function::powerlaw(double x0, double x1, double alpha_g)
{
    //x0 upper,  x1 lower boundary
    if (x0 == 0)
    {
        throw Error_throw_pause("Error! Lower boundary cannot be zero!\nIn class 'Random_function', function 'powerlaw'!\n");
    }
    double y_g = rn(); // random double numbers between 0 and 1
    double x_g = (pow(x1, 1 - alpha_g) - pow(x0, 1 - alpha_g)) * y_g + pow(x0, 1 - alpha_g);
    x_g = pow(x_g, 1 / (1 - alpha_g));
    //cout << x_g << endl;
    return x_g;
}
}; // namespace DFN