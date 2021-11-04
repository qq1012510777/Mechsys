#pragma once
#include "../Error_throw/Error_throw.h"
#include "../Geometry_H/GeoOri_to_NormVec.h"
#include "../Geometry_H/NormVec_to_GeoOri.h"
#include "../Geometry_H/Rotation_verts.h"
#include "../Math_WL_H/erfinv.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include <cmath>
#include <ctime>
#include <iostream>

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
    Random_function(gsl_rng *seed_factor);                           ///< change the "seed_factor" to generate different data series
                                                                     ///< " 0.0 < seed_factor < 1.0"
    double rn();                                                     ///< 0-1 uniform random data
    double unifrm(double a, double b);                               ///< a-b uniform random data
    double gauss(double mean, double sigma, double min, double max); ///< normal distributed random data
    double erlang(double mean, double min, double max);              ///< exponential distributed data

    double powerlaw(double x0, double x1, double alpha_g); ///< powerlaw

    void Fisher_(double mean_y, double mean_x, double fisher_k, ///< Fisher, 1st parameter: dip direction, 2nd: dip angle, 3rd: fisher constant,
                 double min_dip_dir, double max_dip_dir,
                 double min_dip_ang, double max_dip_ang,
                 double &dd_alpha1, double &da_beta1); // use this one that is more accurate!!!

    void Fisher_1(double mean_y, double mean_x, double fisher_k, ///< Fisher, 1st parameter: dip direction, 2nd: dip angle, 3rd: fisher constant,
                  double min_dip_dir, double max_dip_dir,
                  double min_dip_ang, double max_dip_ang,
                  double &dd_alpha1, double &da_beta1); // use this one that is more accurate!!!
    int Bernoulli(double P);

    double lognor_truncated(double mean, double sigma, double min, double max);

    ~Random_function();
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

inline double Random_function::lognor_truncated(double mean, double sigma, double min, double max)
{
    double mu_1 = log(mean * mean / (pow(sigma + mean * mean, 0.5)));
    double sigma_1 = pow(log(1 + ((double)sigma) / (mean * mean)), 0.5); //sigma_1 is input std. deviation

    double l = (log(min) - mu_1) / sigma_1;
    double u = (log(max) - mu_1) / sigma_1;

    double p_l = gsl_cdf_gaussian_P(l, 1);
    double p_u = gsl_cdf_gaussian_P(u, 1);

    double x = rn();

    double x_prime = p_l + (p_u - p_l) * x;

    double z = gsl_cdf_gaussian_Pinv(x_prime, 1);

    double p = exp(z * sigma_1 + mu_1);
    //cout << p << endl;
    return p;
};

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

inline void Random_function::Fisher_(double mean_y, double mean_x, double fisher_k, ///< Fisher, 1st parameter: dip direction, 2nd: dip angle, 3rd: fisher constant,
                                     double min_dip_dir, double max_dip_dir,
                                     double min_dip_ang, double max_dip_ang,
                                     double &dd_alpha1, double &da_beta1)
{
    if (fisher_k < 0 || abs(fisher_k) < 1e-4)
    {
        string AS = "Error! kappa must be larger than zero!\n";
        throw Error_throw_pause(AS);
    }

    size_t n = 0;
    for (size_t i = 0;; ++i)
    {
        n++;
        //--------------------generated orientation along z axis
        double theta = 0, phi = 0;

        phi = this->unifrm(0.0001, 359.9999);
        theta = acos(log(exp(fisher_k) - (exp(fisher_k) - exp(-1.0 * fisher_k)) * this->rn()) / fisher_k);
        theta = theta * 180.0 / M_PI;

        double r = 1;

        double l = r * sin(theta * M_PI / 180.0) * cos(phi * M_PI / 180.0);
        double m = r * sin(theta * M_PI / 180.0) * sin(phi * M_PI / 180.0);
        double n = r * cos(theta * M_PI / 180.0);

        if (n < 0)
        {
            l = -l;
            m = -m;
            n = -n;
        }

        r = (Vector3d{l, m, n}).norm();
        double beta_tmp = acos(n / r) * 180.0 / M_PI;
        double alpha_tmp = atan2(m, l) * 180.0 / M_PI;

        //------now I do not rotate the direction
        if (alpha_tmp < 0)
            alpha_tmp = 360 + alpha_tmp;

        da_beta1 = beta_tmp;

        if (alpha_tmp <= 90)
            dd_alpha1 = 90 - alpha_tmp;

        else if (alpha_tmp > 90)
            dd_alpha1 = 450 - alpha_tmp;

        break;

        /*
        //cout << phi << ", " << theta << endl;
        Vector3d a; // the normal vector of the generated orientation

        a(0) = sin(theta * M_PI / 180) * cos(phi / 180.0 * M_PI);
        a(1) = sin(theta / 180.0 * M_PI) * sin(phi / 180.0 * M_PI);
        a(2) = cos(theta / 180.0 * M_PI);

        if (a(2) < 0)
        {
            a = -a;
        }
        //-------------------------------------------------------

        //----------the mean orientation--------------
        Vector3d b, axis_r;
        double dd = mean_y, da = mean_x;
        Find_normal_vec(dd, da, b);
        //cout << "mean ori: " << b.transpose() << endl;
        //---------------------------------------

        DFN::Vector_2 v{b, axis_r};
        //cout << "axis in xy plane: " << axis_r.transpose() << endl;

        double R_tmp = mean_x * M_PI / 180;
        Quaternion_t Q_axis_1;

        std::vector<Vector3d> f1(1), Verts_1(1);
        f1[0] = a;

        DFN::Rotation_verts R1(f1, -1.0 * R_tmp, Q_axis_1, axis_r, Verts_1);
        // rotate the generated orientation to the mean one

        double l = Verts_1[0][0], m = Verts_1[0][1], n = Verts_1[0][2];
        double norm_ = Verts_1[0].norm();

        if (norm_ < 1e-4)
        {
            string AS = "Error! normal of a frac is (0, 0, 0)!\n";
            throw Error_throw_pause(AS);
        }

        if (Verts_1[0][2] < 0)
        {
            l = -Verts_1[0][0], m = -Verts_1[0][1], n = -Verts_1[0][2];
        }

        double beta_tmp = acos(n) * 180.0 / M_PI;
        double alpha_tmp = atan2(m, l) * 180.0 / M_PI;
        //cout << "after rotate: " << alpha_tmp << ", " << beta_tmp << endl << endl;

        if (alpha_tmp < 0)
            alpha_tmp = 360 + alpha_tmp;

        da_beta1 = beta_tmp;

        if (alpha_tmp <= 90)
            dd_alpha1 = 90 - alpha_tmp;

        else if (alpha_tmp > 90)
            dd_alpha1 = 450 - alpha_tmp;

        if (dd_alpha1 > min_dip_dir && dd_alpha1 < max_dip_dir &&
            da_beta1 > min_dip_ang && da_beta1 < max_dip_ang)
        {
            break;
        }
        */
    }
    //cout << "n = " << n << endl;
};

inline void Random_function::Fisher_1(double mean_y, double mean_x, double fisher_k, ///< Fisher, 1st parameter: dip direction, 2nd: dip angle, 3rd: fisher constant,
                                      double min_dip_dir, double max_dip_dir,
                                      double min_dip_ang, double max_dip_ang,
                                      double &dd_alpha1, double &da_beta1)
{
    if (fisher_k < 0 || abs(fisher_k) < 1e-4)
    {
        string AS = "Error! kappa must be larger than zero!\n";
        throw Error_throw_pause(AS);
    }

    /*
    Vector3d Nor_1, Nor_2, Nor_3, Nor_4, Nor_5;
    DFN::GeoOri_to_NormVec G1{min_dip_dir, min_dip_ang, Nor_1};
    DFN::GeoOri_to_NormVec G2{max_dip_dir, max_dip_ang, Nor_2};
    DFN::GeoOri_to_NormVec G3{min_dip_dir, max_dip_ang, Nor_3};
    DFN::GeoOri_to_NormVec G4{max_dip_dir, min_dip_ang, Nor_4};
    DFN::GeoOri_to_NormVec G5{mean_y, mean_x, Nor_5};

    // rotate normal vector
    Vector3d axis_r1;
    DFN::Vector_2 v1{Nor_5, axis_r1};

    double R_angle1 = mean_x * M_PI / 180;
    Quaternion_t Q_axis;

    std::vector<Vector3d> F1(1), F2(1), F3(1), F4(1),
        E1(1), E2(1), E3(1), E4(1), F5(1), E5(1);

    F1[0] = Nor_1;
    F2[0] = Nor_2;
    F3[0] = Nor_3;
    F4[0] = Nor_4;
    F5[0] = Nor_5;

    DFN::Rotation_verts K1(F1, R_angle1, Q_axis, axis_r1, E1);
    DFN::Rotation_verts K2(F2, R_angle1, Q_axis, axis_r1, E2);
    DFN::Rotation_verts K3(F3, R_angle1, Q_axis, axis_r1, E3);
    DFN::Rotation_verts K4(F4, R_angle1, Q_axis, axis_r1, E4);
    DFN::Rotation_verts K5(F5, R_angle1, Q_axis, axis_r1, E5);

    double dd1 = 0, da1 = 0,
           dd2 = 0, da2 = 0,
           dd3 = 0, da3 = 0,
           dd4 = 0, da4 = 0;

    DFN::NormVec_to_GeoOri P1{E1[0], dd1, da1};
    DFN::NormVec_to_GeoOri P2{E2[0], dd2, da2};
    DFN::NormVec_to_GeoOri P3{E3[0], dd3, da3};
    DFN::NormVec_to_GeoOri P4{E4[0], dd4, da4};

    cout << "1: " << dd1 << ", " << da1 << endl;
    cout << "2: " << dd2 << ", " << da2 << endl;
    cout << "3: " << dd3 << ", " << da3 << endl;
    cout << "4: " << dd4 << ", " << da4 << endl;
    //cout << E5[0].transpose() << endl;
    Vector2d ori1, ori2, ori3, ori4;
    ori1 << dd1, da1;
    ori2 << dd2, da2;
    ori3 << dd3, da3;
    ori4 << dd4, da4;

    if (abs(E5[0][0]) > 1e-4 || abs(E5[0][1]) > 1e-4)
    {
        string AS = "Wrong ratate!\n";
        throw Error_throw_pause(AS);
    }
    */
    size_t n = 0;
    for (size_t i = 0; i < 1; ++i)
    {
        n++;
        //--------------------generated orientation along z axis
        double theta = 0, phi = 0;

        phi = this->unifrm(0.0001, 359.9999);
        //cout << "kappa: " << fisher_k << endl;
        //cout << "upper of theta: " << max_dip_ang << ", lower of theta" << min_dip_ang << endl;

        double C = -1.0 * (exp(fisher_k) - exp(-fisher_k)) / (exp(fisher_k * cos(max_dip_ang * M_PI / 180.0)) - exp(fisher_k * cos(min_dip_ang * M_PI / 180.0)));
        double y_ = this->rn();
        theta = acos((log(-1.0 * y_ * (exp(fisher_k) - exp(-fisher_k)) / C + exp(fisher_k * cos(min_dip_ang * M_PI / 180.0)))) / fisher_k);

        //theta = acos(log(exp(fisher_k) - (exp(fisher_k) - exp(-1.0 * fisher_k)) * this->rn()) / fisher_k);
        theta = theta * 180.0 / M_PI;
        //cout << C << endl;
        //cout << phi << ", " << theta << ", " << y_ << endl;

        if (phi < 0)
            phi = 360 + phi;

        da_beta1 = theta;
        //cout << da_beta1 << endl;
        if (phi <= 90)
            dd_alpha1 = 90 - phi;
        else if (phi > 90)
            dd_alpha1 = 450 - phi;
        /*
        Vector3d a; // the normal vector of the generated orientation

        a(0) = sin(theta * M_PI / 180) * cos(phi / 180.0 * M_PI);
        a(1) = sin(theta / 180.0 * M_PI) * sin(phi / 180.0 * M_PI);
        a(2) = cos(theta / 180.0 * M_PI);

        if (a(2) < 0)
        {
            a = -a;
        }
        //-------------------------------------------------------

        //----------the mean orientation--------------
        Vector3d b, axis_r;
        double dd = mean_y, da = mean_x;
        Find_normal_vec(dd, da, b);
        //cout << "mean ori: " << b.transpose() << endl;
        //---------------------------------------

        DFN::Vector_2 v{b, axis_r};
        //cout << "axis in xy plane: " << axis_r.transpose() << endl;

        double R_tmp = mean_x * M_PI / 180;
        Quaternion_t Q_axis_1;

        std::vector<Vector3d> f1(1), Verts_1(1);
        f1[0] = a;

        DFN::Rotation_verts R1(f1, -1.0 * R_tmp, Q_axis_1, axis_r, Verts_1);
        // rotate the generated orientation to the mean one

        double l = Verts_1[0][0], m = Verts_1[0][1], n = Verts_1[0][2];
        double norm_ = Verts_1[0].norm();

        if (norm_ < 1e-4)
        {
            string AS = "Error! normal of a frac is (0, 0, 0)!\n";
            throw Error_throw_pause(AS);
        }

        if (Verts_1[0][2] < 0)
        {
            l = -Verts_1[0][0], m = -Verts_1[0][1], n = -Verts_1[0][2];
        }

        double beta_tmp = acos(n) * 180.0 / M_PI;
        double alpha_tmp = atan2(m, l) * 180.0 / M_PI;
        //cout << "after rotate: " << alpha_tmp << ", " << beta_tmp << endl << endl;

        if (alpha_tmp < 0)
            alpha_tmp = 360 + alpha_tmp;

        da_beta1 = beta_tmp;

        if (alpha_tmp <= 90)
            dd_alpha1 = 90 - alpha_tmp;

        else if (alpha_tmp > 90)
            dd_alpha1 = 450 - alpha_tmp;

        if (dd_alpha1 > min_dip_dir && dd_alpha1 < max_dip_dir &&
            da_beta1 > min_dip_ang && da_beta1 < max_dip_ang)
        {
            break;
        }*/
    }
    //cout << "n = " << n << endl;
};

inline int Random_function::Bernoulli(double P)
{
    //cout << gsl_ran_bernoulli(dmodul, P) << endl;
    return gsl_ran_bernoulli(dmodul, P);
};

inline Random_function::~Random_function(){
    //gsl_rng_free(dmodul);
};

}; // namespace DFN