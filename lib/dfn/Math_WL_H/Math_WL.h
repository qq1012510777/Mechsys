#pragma once
#include "../Eigen_API/Eigen_API.h"
#include "../Line_WL_H/Line_WL.h"
#include "../Quaternion_H/Quaternion.h"
#include "Eigen/Dense"
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <vector>
using namespace std;
using namespace Eigen;

typedef std::chrono::time_point<std::chrono::_V2::steady_clock, std::chrono::duration<long int, std::ratio<1, 1000000000>>> Time_Pnt;
void time_counter_start(Time_Pnt &start);
void time_counter_end(Time_Pnt start,
                      Time_Pnt &end,
                      string content,
                      string unit);

string To_string_with_width(size_t val, size_t width);

template <typename Ts>
std::string to_string_with_precision(const Ts a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
};

// LAPACK
extern "C"
{
    // DGETRF - compute an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges
    void dgetrf_(const int *m, const int *n, double *a, const int *lda, int *ipiv, int *info);

    // DGETRI - compute the inverse of a matrix using the LU factorization computed by DGETRF
    void dgetri_(const int *n, double *a, const int *lda, int *ipiv, double *work, const int *lwork, int *info);

    // DGESVD - computes the singular value decomposition of a real M-by-N matrix A, optionally computing the left and/or right singular vectors
    void dgesvd_(const char *jobu, const char *jobvt, const int *M, const int *N, double *A, const int *lda, double *S, double *U, const int *ldu, double *VT, const int *ldvt, double *work, const int *lwork, const int *info);

    // DGESV - double-general-solver
    void dgesv_(int *Np, int *NRHSp, double *A, int *LDAp, int *IPIVp, double *B, int *LDBp, int *INFOp);

    // DGEMV - double-matrix-vector
    void dgemv_(const char *jobu, const int *m, const int *n, const double *a, const double *A, const int *lda, const double *x, const int *incx, const double *b, double *y, const int *incy);
}

//< declaration

/// based on orientation, find normal vector
void Find_normal_vec(double dip_direction,
                     double dip_angle,
                     Vector3d &a);

// used in function: random_unsigned_integer
#define Random_1(low, up) (rand() % (up - low + 1)) + low

size_t random_unsigned_integer(size_t x, size_t y);

//generate random numbers
double random_double(size_t min, size_t max);

//the function below finds the vector that
//(1) is vertical to fracture normal vector;
//and (2) lies on the horizontal plane (z = 0)
void Find_vector_2(Vector3d Normal_vector,
                   Vector3d &temp3);

//the function below
//determine whether two
//infinite planes are parallel or not
//when index1 = 1, means parallel;
//if index1 = 0, means two infinite planes intersect
//if index2 = 0 means parallel but not overlap,
//if index2 = 1 means overlap
void Parallel_or_not(Vector4d plane_parameter1,
                     Vector4d plane_parameter2,
                     size_t &index1,
                     size_t &index2);

//the function below finds the maximum element
double Find_max_z_value(std::vector<Vector3d> A);

//below finds minimum element
double Find_min_z_value(std::vector<Vector3d> A);

//below finds the intersection between
//a fracture and a infinite plane
//(this plane must be horizontal)
void Intersection_between_2ndpolygon_and_infinite_plane(size_t vertexes_of_polygon_,
                                                        double z_zplane,
                                                        std::vector<Vector3d> array,
                                                        std::vector<Vector3d> &Intersection_1);

//below finds a relatively infinte line, i
//.e., absolute values of coordinates of
//both two ends are very large
//also they are in horizontal plane,
//or say, 2D Cartesian system
void Output_a_relatively_infinite_line(double max,
                                       std::vector<Vector3d> Intersection_1,
                                       std::vector<Vector3d> &Intersection_infinite);

//below finds intersection between line segment
//and a polygon, also, in 2D Cartesian system
///there are in the same horizontal plane
void Intersection_between_line_segment_and_polygon(size_t &numOfIntersectionPoint,
                                                   std::vector<Vector3d> Verts_1,
                                                   std::vector<Vector3d> Intersection_infinite,
                                                   std::vector<Vector3d> &Intersection_2);

//below determines whether two lines
//are intersect, also, in 2D Cartesian system
//determine if two lines intersect
size_t is_intersect(DFN::Line myline1,
                    DFN::Line myline2);

//below finds the intersection between two line segments
//(if they have)
void intersection_between_two_line_segments(std::vector<Vector3d> Intersection_infinite,
                                            Vector3d A,
                                            Vector3d B,
                                            double &intersection_x,
                                            double &intersection_y);

void intersection_between_two_line_segments(const double x0,
                                            const double y0,
                                            const double x1,
                                            const double y1,
                                            const double x2,
                                            const double y2,
                                            const double x3,
                                            const double y3,
                                            double &intersection_x,
                                            double &intersection_y);

//below finds the intersection between
//two 1D intervals (if they have)
void Intersection_of_1D_intervals(std::vector<Vector3d> Intersection_4,
                                  std::vector<Vector3d> Intersection_5,
                                  std::vector<Vector3d> &Intersection_6);

/*****************************************************************
the function below can be used to calculate first Modified Bessel function (when the input, i.e., n, is 0, 1, 2, 3, or 4)
they are equal to I0(x), I1(x), ......
******************************************************************/
double first_modified_Bessel(int n, double x);

typedef long long LL;
double round(double number, unsigned int bits);
// round a number

bool Find_repetitive_thing(const Vector3d A,
                           const std::vector<std::vector<Vector3d>> JXY,
                           const size_t x,
                           const size_t y,
                           int &x1,
                           int &y1);

bool Find_repetitive_thing_p(const Vector3d A,
                             const std::vector<std::vector<Vector3d>> JXY,
                             const size_t x,
                             const size_t y,
                             int &x1,
                             int &y1,
                             std::vector<size_t> NO_Nodes_p_each_frac);

//

void binary_linear_equation(const Vector3d A,
                            const Vector3d B,
                            const int N /*NO of subsegments*/,
                            std::vector<Vector3d> &H);

void binary_linear_equation_Trans(const Vector3d A,
                                  const Vector3d B,
                                  Vector3d &C,
                                  const double xi,
                                  const double X1 = -1,
                                  const double X2 = 1);

bool comp(std::pair<size_t, double> &a, std::pair<size_t, double> &b);

void Sort_trace_intersection_pnts(std::vector<Vector3d> &Intersec_traces_sub,
                                  const Vector3d A);
// sort, the 2nd argument is the initial pnt

double Length_2d_segment(const Vector3d A, const Vector3d B);

double unkown_vector_operator(const Vector3d p1, const Vector3d p2);
//

bool pointInConvexPolygon(const Vector3d p, const std::vector<Vector3d> polygon);

// the PNT must be inside the domain
bool Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(const Vector3d A, const Vector3d B, const Vector3d C, const Vector3d D, const Vector3d PNT);

bool Determine_IF_a_pnt_lies_in_a_2D_line_segment(const Vector3d A, const Vector3d B, const Vector3d C);

double Numerical_integration_linear_3(VectorXd Vec_q_normal, double L);

double Numerical_integration_linear_2(VectorXd Vec_q_normal, double L);

double Heron_formula(Vector3d pnt0, Vector3d pnt2, Vector3d pnt4);

void intersection_between_two_overlaped_line_segments(std::vector<Vector3d> Intersection_infinite,
                                                      Vector3d A,
                                                      Vector3d B,
                                                      std::vector<Vector3d> &Intersection_2);

//*****************************************

//< definition

void Find_normal_vec(double dip_direction,
                     double dip_angle,
                     Vector3d &a)
{
    ///spherical system firstly
    double alpha = 0, beta = 0;
    beta = dip_angle;
    if (dip_direction >= 90)
        alpha = 450 - dip_direction;
    else if (dip_direction <= 90)
        alpha = 90 - dip_direction;

    //------------------
    a(0) = sin(beta * M_PI / 180) * cos(alpha / 180.0 * M_PI);
    a(1) = sin(beta / 180.0 * M_PI) * sin(alpha / 180.0 * M_PI);
    a(2) = cos(beta / 180.0 * M_PI);

    if (a(2) < 0)
    {
        a = -a;
    }
};

size_t random_unsigned_integer(size_t x, size_t y)
{
    return Random_1(x, y);
};

double random_double(size_t min, size_t max) //generate random numbers
{
    double m1 = (double)(rand() % 101) / 101;
    min++;
    double m2 = (double)((rand() % (max - min + 1)) + min);
    m2 = m2 - 1;
    return m1 + m2;
};

void Find_vector_2(Vector3d Normal_vector, Vector3d &temp3)
{
    temp3(2) = 0;
    if (Normal_vector(0) > 0)
    {
        if (Normal_vector(1) > 0)
        {
            temp3(0) = -Normal_vector(1);
            temp3(1) = Normal_vector(0);
        }
        else if (Normal_vector(1) < 0)
        {
            temp3(0) = Normal_vector(1);
            temp3(1) = -Normal_vector(0);
        }
        else if (Normal_vector(1) == 0)
        {
            temp3(0) = 0;
            temp3(1) = Normal_vector(0);
        }
    }
    else if (Normal_vector(0) < 0)
    {
        if (Normal_vector(1) < 0)
        {
            temp3(0) = -Normal_vector(1);
            temp3(1) = Normal_vector(0);
        }
        else if (Normal_vector(1) > 0)
        {
            temp3(0) = -Normal_vector(1);
            temp3(1) = Normal_vector(0);
        }
        else if (Normal_vector(1) == 0)
        {
            temp3(0) = 0;
            temp3(1) = Normal_vector(0);
        }
    }
    else if (Normal_vector(0) == 0)
    {
        if (Normal_vector(1) > 0)
        {
            temp3(0) = -Normal_vector(1);
            temp3(1) = 0;
        }
        else if (Normal_vector(1) < 0)
        {
            temp3(0) = -Normal_vector(1);
            temp3(1) = 0;
        }
        else if (Normal_vector(1) == 0)
        {
            temp3(0) = 0;
            temp3(1) = 0;
            return;
        }
    };
}

void Parallel_or_not(Vector4d plane_parameter1,
                     Vector4d plane_parameter2,
                     size_t &index1,
                     size_t &index2)
{
    double a1 = plane_parameter1(0);
    double b1 = plane_parameter1(1);
    double c1 = plane_parameter1(2);
    double d1 = plane_parameter1(3);
    double a2 = plane_parameter2(0);
    double b2 = plane_parameter2(1);
    double c2 = plane_parameter2(2);
    double d2 = plane_parameter2(3);

    index1 = 0; //when index1 = 1, means parallel; if index1 = 0, means two infinite planes intersect
    index2 = 0; //0 means parallel but not overlap, 1 means overlap

    Vector3d A, B;
    A << a1, b1, c1;
    B << a2, b2, c2;

    if (abs((A.cross(B)).norm()) < 0.0001)
    {
        index1 = 1;
    }
    else
    {
        index1 = 0;
    }

    if (index1 == 1)
    {
        double x = 1.5, y = 3.6, z1 = 0, z2 = 0;
        z1 = -(a1 * x + b1 * y + d1) / c1;
        z2 = -(a2 * x + b2 * y + d2) / c2;
        if (abs(z1 - z2) < 0.001)
        {
            index2 = 1;
        }
        else
            index2 = 0;
    }
};

//the function below finds the maximum element
double Find_max_z_value(std::vector<Vector3d> A)
{

    double Max = 0;
    for (size_t i = 0; i < A.size() - 1; i++)
    {
        if (i == 0)
            Max = max(A[0](2), A[1](2));
        else
            Max = max(Max, A[i + 1](2));
    };
    return Max;
};

//below finds minimum element
double Find_min_z_value(std::vector<Vector3d> A)
{

    double Min = 0;
    for (size_t i = 0; i < A.size() - 1; i++)
    {
        if (i == 0)
            Min = min(A[0](2), A[1](2));
        else
            Min = min(Min, A[i + 1](2));
    }
    return Min;
};

void Intersection_between_2ndpolygon_and_infinite_plane(size_t vertexes_of_polygon_,
                                                        double z_zplane,
                                                        std::vector<Vector3d> array,
                                                        std::vector<Vector3d> &Intersection_1)
{

    double end_point1_x = 0;
    double end_point1_y = 0;
    double end_point2_x = 0;
    double end_point2_y = 0;

    size_t k = 0;
    for (size_t i = 0; i < vertexes_of_polygon_; i++)
    {
        if (i == vertexes_of_polygon_ - 1)
        {
            if ((array[i](2) <= z_zplane && z_zplane <= array[i + 1](2)) || (array[i](2) >= z_zplane && z_zplane >= array[i + 1](2)))
            {
                double t = (z_zplane - array[0](2)) / (array[i](2) - array[0](2));
                end_point1_x = array[0](0) + (array[i](0) - array[0](0)) * t;
                end_point1_y = array[0](1) + (array[i](1) - array[0](1)) * t;
                k = i + 1;

                break;
            }
        }
        else
        {
            if ((array[i](2) <= z_zplane && z_zplane <= array[i + 1](2)) || (array[i](2) >= z_zplane && z_zplane >= array[i + 1](2)))
            {
                double t = (z_zplane - array[i + 1](2)) / (array[i](2) - array[i + 1](2));
                end_point1_x = array[i + 1](0) + (array[i](0) - array[i + 1](0)) * t;
                end_point1_y = array[i + 1](1) + (array[i](1) - array[i + 1](1)) * t;
                k = i + 1;

                break;
            }
        }
    }

    for (size_t j = k; j < vertexes_of_polygon_; j++)
    {

        if (j == vertexes_of_polygon_ - 1)
        {

            if ((array[j](2) <= z_zplane && z_zplane <= array[0](2)) || (array[j](2) >= z_zplane && z_zplane >= array[0](2)))
            {
                double t = (z_zplane - array[0](2)) / (array[j][2] - array[0](2));

                end_point2_x = array[0](0) + (array[j](0) - array[0](0)) * t;
                end_point2_y = array[0](1) + (array[j](1) - array[0](1)) * t;

                j = vertexes_of_polygon_;
            }
        }
        else
        {
            if ((array[j](2) <= z_zplane && z_zplane <= array[j + 1](2)) || (array[j](2) >= z_zplane && z_zplane >= array[j + 1](2)))
            {
                double t = (z_zplane - array[j + 1](2)) / (array[j](2) - array[j + 1](2));

                end_point2_x = array[j + 1](0) + (array[j](0) - array[j + 1](0)) * t;
                end_point2_y = array[j + 1](1) + (array[j](1) - array[j + 1](1)) * t;

                j = vertexes_of_polygon_;
            }
        }
        if (j == vertexes_of_polygon_ - 1)
            break;
    }
    Intersection_1.resize(2);
    Intersection_1[0] << end_point1_x, end_point1_y, 0;
    Intersection_1[1] << end_point2_x, end_point2_y, 0;
}

void Output_a_relatively_infinite_line(double max,
                                       std::vector<Vector3d> Intersection_1,
                                       std::vector<Vector3d> &Intersection_infinite)
{
    double x0, y0, x1, y1, x2, y2, x3, y3;
    x0 = Intersection_1[0](0);
    y0 = Intersection_1[0](1);
    x1 = Intersection_1[1](0);
    y1 = Intersection_1[1](1);
    x2 = max;
    x3 = -max;
    if ((x1 - x0) != 0)
    {
        double t1 = (x2 - x0) / (x1 - x0);
        y2 = t1 * (y1 - y0) + y0;

        double t2 = (x3 - x0) / (x1 - x0);
        y3 = t2 * (y1 - y0) + y0;
    }
    else
    {
        y2 = x2;
        y3 = x3;
        x2 = x0;
        x3 = x0;
    }
    Intersection_infinite.resize(2);
    Intersection_infinite[0](0) = x2;
    Intersection_infinite[0](1) = y2;
    Intersection_infinite[1](0) = x3;
    Intersection_infinite[1](1) = y3;
};

size_t is_intersect(DFN::Line myline1,
                    DFN::Line myline2)
{
    double x0, y0, x1, y1, x2, y2, x3, y3;
    x0 = myline1.xa;
    y0 = myline1.ya;
    x1 = myline1.xb;
    y1 = myline1.yb;
    x2 = myline2.xa;
    y2 = myline2.ya;
    x3 = myline2.xb;
    y3 = myline2.yb;

    //move to origin
    x1 -= x0;
    y1 -= y0;
    x2 -= x0;
    y2 -= y0;
    x3 -= x0;
    y3 -= y0;
    x0 = 0;
    y0 = 0;

    //rotation to x-axis
    Vector3d axis_z;
    axis_z << 0, 0, -1;
    double R_angle_temp1 = atan2(y1, x1);

    Quaternion_t Q_axis_1;
    NormalizeRotation(R_angle_temp1, axis_z, Q_axis_1);

    Vector3d temp1;
    temp1 << x1, y1, 0;
    Vector3d temp2;
    Rotation(temp1, Q_axis_1, temp2);
    x1 = temp2(0);
    y1 = temp2(1);

    Vector3d temp3, temp4;
    temp3 << x2, y2, 0;
    Rotation(temp3, Q_axis_1, temp4);
    x2 = temp4(0);
    y2 = temp4(1);

    Vector3d temp5, temp6;
    temp5 << x3, y3, 0;
    Rotation(temp5, Q_axis_1, temp6);
    x3 = temp6(0);
    y3 = temp6(1);

    Vector3d A1, A2, B1, B2;
    A1 << x0, y0, 0;
    A2 << x1, y1, 0;

    B1 << x2, y2, 0;
    B2 << x3, y3, 0;
    DFN::Line A(A1, A2), B(B1, B2);

    if ((abs(y0) < 0.0001 && abs(y1) < 0.0001 && abs(y2) < 0.0001 && abs(y3) < 0.0001) &&
        ((A.get_max_x() - B.get_min_x() > -0.0001) || (B.get_max_x() - A.get_min_x() > -0.0001)))
    {
        //cout << "overlapping;\n";
        return 3;
    }

    if ((A.get_max_x() - B.get_min_x() < -0.0001) ||
        (B.get_max_x() - A.get_min_x() < -0.0001) ||
        (A.get_max_y() - B.get_min_y() < -0.0001) ||
        (B.get_max_y() - A.get_min_y() < -0.0001))
    {
        //cout << " 000 disconnected;\n";
        return 0;
    } //0 means disconnected
    double res1 = (A.xa - A.xb) * (B.ya - A.yb) - (A.ya - A.yb) * (B.xa - A.xb);
    double res2 = (A.xa - A.xb) * (B.yb - A.yb) - (A.ya - A.yb) * (B.xb - A.xb);

    double res3 = (B.xa - B.xb) * (A.ya - B.yb) - (B.ya - B.yb) * (A.xa - B.xb);
    double res4 = (B.xa - B.xb) * (A.yb - B.yb) - (B.ya - B.yb) * (A.xb - B.xb);
    if (res1 * res2 <= 0.0001 && res3 * res4 <= 0.0001)
    {
        //cout << " 111 connected;\n";
        //cout << res1 * res2 << ", " << res3 * res4 << endl;

        //cout << "intersected;\n";
        return 1;

    } //1 means connected
    else
    {
        //cout << " 111 disconnected;\n";
        //cout << res1 * res2 << ", " << res3 * res4 << endl;
        return 0;
    }
};

void intersection_between_two_line_segments(std::vector<Vector3d> Intersection_infinite,
                                            Vector3d A,
                                            Vector3d B,
                                            double &intersection_x,
                                            double &intersection_y)
{
    //// if two lines are overlapping, how to do??????????
    double x0 = Intersection_infinite[0](0);
    double y0 = Intersection_infinite[0](1);
    double x1 = Intersection_infinite[1](0);
    double y1 = Intersection_infinite[1](1);
    double x2 = A(0);
    double y2 = A(1);
    double x3 = B(0);
    double y3 = B(1);

    if (abs(x1 - x0) > 0.001 && abs(x3 - x2) > 0.001)
    {
        double k1, k2;
        k1 = round((y1 - y0), 4) / (x1 - x0);
        k2 = round((y3 - y2), 4) / (x3 - x2);
        if (k1 != 0 && k2 != 0)
        {
            intersection_x = (y0 - y2 + k2 * x2 - k1 * x0) / (k2 - k1);
            intersection_y = k1 * (intersection_x - x0) + y0;
        }
        else if ((k1 == 0 && k2 != 0) || (k1 != 0 && k2 == 0))
        {
            if (k1 == 0)
            {
                intersection_y = y0;
                intersection_x = (intersection_y - y2) / ((y3 - y2) / (x3 - x2)) + x2;
            }
            else if (k2 == 0)
            {
                intersection_y = y2;
                intersection_x = (intersection_y - y0) / ((y1 - y0) / (x1 - x0)) + x0;
            }
        }
        else
        {
            //two overlapping and it should be excluded by is_intersect function
            cout << "two horizontal lines overlap, the function 'is_intersect' should be used beforehand!\n";
            exit(0);
        }
    }
    else if ((abs(x1 - x0) < 0.001 && abs(x3 - x2) > 0.001) || (abs(x1 - x0) > 0.001 && abs(x3 - x2) < 0.001))
    {
        if (abs(x1 - x0) < 0.001 && abs(x3 - x2) > 0.001)
        {
            intersection_x = x0;
            double k2 = (y3 - y2) / (x3 - x2);
            intersection_y = k2 * (intersection_x - x2) + y2;
        }
        else if (abs(x1 - x0) > 0.001 && abs(x3 - x2) < 0.001)
        {
            intersection_x = x2;
            double k1 = (y1 - y0) / (x1 - x0);
            intersection_y = k1 * (intersection_x - x0) + y0;
        }
    }
    else
    {
        //two overlapping and it should be excluded by is_intersect function
        cout << "two vertical lines overlap, the function 'is_intersect' should be used beforehand!\n";
        exit(0);
    }
}

void intersection_between_two_line_segments(const double x0,
                                            const double y0,
                                            const double x1,
                                            const double y1,
                                            const double x2,
                                            const double y2,
                                            const double x3,
                                            const double y3,
                                            double &intersection_x,
                                            double &intersection_y)
{

    if (abs(x1 - x0) > 0.001 && abs(x3 - x2) > 0.001)
    {
        double k1, k2;
        k1 = round((y1 - y0), 4) / (x1 - x0);
        k2 = round((y3 - y2), 4) / (x3 - x2);
        if (k1 != 0 && k2 != 0)
        {
            intersection_x = (y0 - y2 + k2 * x2 - k1 * x0) / (k2 - k1);
            intersection_y = k1 * (intersection_x - x0) + y0;
        }
        else if ((k1 == 0 && k2 != 0) || (k1 != 0 && k2 == 0))
        {
            if (k1 == 0)
            {
                intersection_y = y0;
                intersection_x = (intersection_y - y2) / ((y3 - y2) / (x3 - x2)) + x2;
            }
            else if (k2 == 0)
            {
                intersection_y = y2;
                intersection_x = (intersection_y - y0) / ((y1 - y0) / (x1 - x0)) + x0;
            }
        }
        else
        {
            //two overlapping and it should be excluded by is_intersect function
            cout << "two horizontal lines overlap, the function 'is_intersect' should be used beforehand!\n";
            exit(0);
        }
    }
    else if ((abs(x1 - x0) < 0.001 && abs(x3 - x2) > 0.001) || (abs(x1 - x0) > 0.001 && abs(x3 - x2) < 0.001))
    {
        if (abs(x1 - x0) < 0.001 && abs(x3 - x2) > 0.001)
        {
            intersection_x = x0;
            double k2 = (y3 - y2) / (x3 - x2);
            intersection_y = k2 * (intersection_x - x2) + y2;
        }
        else if (abs(x1 - x0) > 0.001 && abs(x3 - x2) < 0.001)
        {
            intersection_x = x2;
            double k1 = (y1 - y0) / (x1 - x0);
            intersection_y = k1 * (intersection_x - x0) + y0;
        }
    }
    else
    {
        //two overlapping and it should be excluded by is_intersect function
        cout << "two vertical lines overlap, the function 'is_intersect' should be used beforehand!\n";
        exit(0);
    }
}

void Intersection_between_line_segment_and_polygon(size_t &numOfIntersectionPoint,
                                                   std::vector<Vector3d> Verts_1,
                                                   std::vector<Vector3d> Intersection_infinite,
                                                   std::vector<Vector3d> &Intersection_2)
{

    double end_point3_x = 0;
    //points between infinite line and
    //each edge of 1st fracture
    double end_point3_y = 0;
    double end_point4_x = 0;
    double end_point4_y = 0;

    size_t index_tmp2 = 0;
    size_t index_tmp3 = 0;
    for (size_t i = 0; i < Verts_1.size(); i++)
    {
        index_tmp3 = i;
        if (i != Verts_1.size() - 1)
        {
            DFN::Line l1(Intersection_infinite);
            DFN::Line l2(Verts_1[i], Verts_1[i + 1]);
            size_t r1 = is_intersect(l1, l2);
            if (r1 == 1)
            {
                index_tmp2++;
                intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[i + 1], end_point3_x, end_point3_y);
                break;
            }
            else if (r1 == 3)
            {
                Intersection_2.resize(2);
                intersection_between_two_overlaped_line_segments(Intersection_infinite, Verts_1[i], Verts_1[i + 1], Intersection_2);
                return;
            }
        }
        else
        {

            DFN::Line l1(Intersection_infinite);
            DFN::Line l2(Verts_1[i], Verts_1[0]);
            size_t r1 = is_intersect(l1, l2);
            if (r1 == 1)
            {
                index_tmp2++;
                intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[0], end_point3_x, end_point3_y);

                break;
            }
            else if (r1 == 3)
            {
                Intersection_2.resize(2);
                intersection_between_two_overlaped_line_segments(Intersection_infinite, Verts_1[i], Verts_1[0], Intersection_2);
                return;
            }
        }
    }
    for (size_t i = index_tmp3 + 1; i < Verts_1.size(); i++)
    {
        if (i != Verts_1.size() - 1)
        {
            DFN::Line l1(Intersection_infinite);
            DFN::Line l2(Verts_1[i], Verts_1[i + 1]);
            size_t r1 = is_intersect(l1, l2);
            if (r1 == 1)
            {
                index_tmp2++;
                intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[i + 1], end_point4_x, end_point4_y);
                break;
            }
            else if (r1 == 3)
            {
                Intersection_2.resize(2);
                intersection_between_two_overlaped_line_segments(Intersection_infinite, Verts_1[i], Verts_1[i + 1], Intersection_2);
                return;
            }
        }
        else
        {
            DFN::Line l1(Intersection_infinite);
            DFN::Line l2(Verts_1[i], Verts_1[0]);
            size_t r1 = is_intersect(l1, l2);

            if (r1 == 1)
            {
                index_tmp2++;
                intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[0], end_point4_x, end_point4_y);
                break;
            }
            else if (r1 == 3)
            {
                Intersection_2.resize(2);
                intersection_between_two_overlaped_line_segments(Intersection_infinite, Verts_1[i], Verts_1[0], Intersection_2);
                return;
            }
        }
    }
    Intersection_2.resize(2);
    Intersection_2[0] << end_point3_x, end_point3_y, 0;
    Intersection_2[1] << end_point4_x, end_point4_y, 0;
    numOfIntersectionPoint = index_tmp2;
}

//below finds the intersection between
//two 1D intervals (if they have)
void Intersection_of_1D_intervals(std::vector<Vector3d> Intersection_4,
                                  std::vector<Vector3d> Intersection_5,
                                  std::vector<Vector3d> &Intersection_6)
{
    //first two are intersection points between horizontal plane and 2nd fracture
    double x1, x2, x3, x4, x5 = 0, x6 = 0;
    x1 = Intersection_4[0](0);
    x2 = Intersection_4[1](0);
    x3 = Intersection_5[0](0);
    x4 = Intersection_5[1](0);

    double low1, up1, low2, up2;
    if (x2 >= x1)
    {
        up1 = x2;
        low1 = x1;
    }
    else
    {
        up1 = x1;
        low1 = x2;
    }
    if (x4 >= x3)
    {
        up2 = x4;
        low2 = x3;
    }
    else
    {
        up2 = x3;
        low2 = x4;
    }
    if (low1 > up2 || low2 > up1)
    {
        //printf("NULL Intersection of two 1D intervals\n");
        return;
    }
    else if (up1 >= up2)
    {
        if (low1 > low2)
        {
            x5 = low1;
            x6 = up2;
        }
        else
        {
            x5 = low2;
            x6 = up2;
        }
    }
    else if (up2 > up1)
    {
        if (low2 > low1)
        {
            x5 = low2;
            x6 = up1;
        }
        else
        {
            x5 = low1;
            x6 = up1;
        }
    }
    Intersection_6[0](0) = x5;
    Intersection_6[0](1) = Intersection_4[0](1);
    Intersection_6[1](0) = x6;
    Intersection_6[1](1) = Intersection_4[0](1);
};

double first_modified_Bessel(int n, double x)
{
    int i, m;
    double t, y, p = 0, b0, b1, q;
    Eigen::VectorXd a(7);
    a << 1.0, 3.5156229, 3.0899424, 1.2067492,
        0.2659732, 0.0360768, 0.0045813;
    Eigen::VectorXd b(7);
    b << 0.5, 0.87890594, 0.51498869,
        0.15084934, 0.02658773, 0.00301532, 0.00032411;
    Eigen::VectorXd c(9);
    c << 0.39894228, 0.01328592, 0.00225319,
        -0.00157565, 0.00916281, -0.02057706,
        0.02635537, -0.01647633, 0.00392377;
    Eigen::VectorXd d(9);
    d << 0.39894228, -0.03988024, -0.00362018,
        0.00163801, -0.01031555, 0.02282967,
        -0.02895312, 0.01787654, -0.00420059;
    if (n < 0)
        n = -n;
    t = fabs(x);
    if (n != 1)
    {
        if (t < 3.75)
        {
            y = (x / 3.75) * (x / 3.75);
            p = a[6];
            for (i = 5; i >= 0; i--)
                p = p * y + a[i];
        }
        else
        {
            y = 3.75 / t;
            p = c[8];
            for (i = 7; i >= 0; i--)
                p = p * y + c[i];
            p = p * exp(t) / sqrt(t);
        }
    }
    if (n == 0)
        return (p);
    q = p;
    if (t < 3.75)
    {
        y = (x / 3.75) * (x / 3.75);
        p = b[6];
        for (i = 5; i >= 0; i--)
            p = p * y + b[i];
        p = p * t;
    }
    else
    {
        y = 3.75 / t;
        p = d[8];
        for (i = 7; i >= 0; i--)
            p = p * y + d[i];
        p = p * exp(t) / sqrt(t);
    }
    if (x < 0.0)
        p = -p;
    if (n == 1)
        return (p);
    if (x == 0.0)
        return (0.0);
    y = 2.0 / t;
    t = 0.0;
    b1 = 1.0;
    b0 = 0.0;
    m = n + (int)sqrt(40.0 * n);
    m = 2 * m;

    for (i = m; i > 0; i--)
    {
        p = b0 + i * y * b1;
        b0 = b1;
        b1 = p;
        if (fabs(b1) > 1.0e+10)
        {
            t = t * 1.0e-10;
            b0 = b0 * 1.0e-10;
            b1 = b1 * 1.0e-10;
        }
        if (i == n)
            t = b0;
    }
    p = t * q / b1;
    if ((x < 0.0) && (n % 2 == 1))
        p = -p;
    return (p);
};

double round(double number, unsigned int bits)
{
    LL integerPart = number;
    number -= integerPart;
    for (unsigned int i = 0; i < bits; ++i)
        number *= 10;
    number = (LL)(number + 0.5);
    for (unsigned int i = 0; i < bits; ++i)
        number /= 10;
    return integerPart + number;
};

bool Find_repetitive_thing(const Vector3d A,
                           const std::vector<std::vector<Vector3d>> JXY,
                           const size_t x,
                           const size_t y,
                           int &x1,
                           int &y1)
{
    x1 = -1;
    y1 = -1;
    for (size_t i = 0; i <= x; ++i)
    {
        for (size_t j = 0; j < JXY[i].size(); ++j)
        {
            if (i == x && j == y)
            {
                return false;
            }
            Vector3d S = A - JXY[i][j];
            if (abs(S(0)) < 0.001 && abs(S(1)) < 0.001 && abs(S(2)) < 0.001)
            {
                x1 = i;
                y1 = j;
                return true;
            }
        }
    }
    return false;
};

bool Find_repetitive_thing_p(const Vector3d A,
                             const std::vector<std::vector<Vector3d>> JXY,
                             const size_t x,
                             const size_t y,
                             int &x1,
                             int &y1,
                             std::vector<size_t> NO_Nodes_p_each_frac)
{
    x1 = -1;
    y1 = -1;
    for (size_t i = 0; i <= x; ++i)
    {
        for (size_t j = 0; j < NO_Nodes_p_each_frac[i]; ++j)
        {
            if (i == x && j == y)
            {
                return false;
            }
            Vector3d S = A - JXY[i][j];
            if (abs(S(0)) < 0.001 && abs(S(1)) < 0.001 && abs(S(2)) < 0.001)
            {
                x1 = i;
                y1 = j;
                return true;
            }
        }
    }
    return false;
};

void binary_linear_equation(const Vector3d A,
                            const Vector3d B,
                            const int N /*NO of subsegments*/,
                            std::vector<Vector3d> &H)
{
    double x0, y0, x1, y1, x2, y2;
    x0 = A(0);
    y0 = A(1);
    x1 = B(0);
    y1 = B(1);

    if (N == 1)
        return;

    if ((abs(x0 - x1) > 0.0001 && abs(y0 - y1) > 0.0001) || (abs(x0 - x1) > 0.0001 && abs(y0 - y1) < 0.0001))
    {
        Vector3d tmp_1, tmp_2;
        tmp_1 << x1 - x0, y1 - y0, 0;
        // cout << "after move A to origin, B = " << tmp_1 << endl;
        double angle_tmp = atan2(y1 - y0, x1 - x0);
        //cout << "angle_tmp = " << 180 / M_PI * angle_tmp << endl;
        Vector3d axis_z;
        axis_z << 0, 0, -1;

        Quaternion_t Q_axis_z;

        NormalizeRotation(angle_tmp, axis_z, Q_axis_z);
        Rotation(tmp_1, Q_axis_z, tmp_2);
        //cout << "after rotation, B = " << tmp_2 << endl;

        double delta_x = (tmp_2(0) - 0) / N;
        for (size_t i = 0; i < (size_t)N - 1; ++i)
        {
            x2 = 0 + (i + 1) * delta_x;
            Quaternion_t Q_axis_z_A;

            Vector3d tmp_A1, tmp_A2;
            tmp_A1 << x2, 0, 0;
            tmp_A2 << 0, 0, 0;
            axis_z << 0, 0, 1;
            double angle_tmp_back = angle_tmp;
            NormalizeRotation(angle_tmp_back, axis_z, Q_axis_z_A);
            Rotation(tmp_A1, Q_axis_z_A, tmp_A2);

            tmp_A2(0) += x0;
            tmp_A2(1) += y0;
            tmp_A2(2) = 0;
            H.push_back(tmp_A2);
        }
    }
    else if (abs(x0 - x1) < 0.0001 && abs(y0 - y1) > 0.0001) // x values are the same
    {
        for (size_t i = 0; i < (size_t)N - 1; ++i)
        {
            double delta_y = (y1 - y0) / N;
            y2 = y0 + (i + 1) * delta_y;
            x2 = x0;
            Vector3d tmp_u;
            tmp_u << x2, y2, 0;
            H.push_back(tmp_u);
        }
    }
    else
    {
        //(abs(x0 - x1) < 0.0001 && abs(y0 - y1) < 0.0001)
        if (H.size() != 0)
            H.resize(0);
        return;
    }
};

void binary_linear_equation_Trans(const Vector3d A,
                                  const Vector3d B,
                                  Vector3d &C,
                                  const double xi,
                                  const double X1,
                                  const double X2)
{
    double x0, y0, x1, y1, x2, y2;
    x0 = A(0);
    y0 = A(1);
    x1 = B(0);
    y1 = B(1);

    if ((abs(x0 - x1) > 0.0001 && abs(y0 - y1) > 0.0001) || (abs(x0 - x1) > 0.0001 && abs(y0 - y1) < 0.0001))
    {
        Vector3d tmp_1, tmp_2;
        tmp_1 << x1 - x0, y1 - y0, 0;
        // cout << "after move A to origin, B = " << tmp_1 << endl;
        double angle_tmp = atan2(y1 - y0, x1 - x0);
        //cout << "angle_tmp = " << 180 / M_PI * angle_tmp << endl;
        Vector3d axis_z;
        axis_z << 0, 0, -1;

        Quaternion_t Q_axis_z;

        NormalizeRotation(angle_tmp, axis_z, Q_axis_z);
        Rotation(tmp_1, Q_axis_z, tmp_2);
        //cout << "after rotation, B = " << tmp_2 << endl;

        double L = tmp_2(0);
        double delta_x = (xi - X1) / (X2 - X1) * L;

        x2 = 0 + (0 + 1) * delta_x;
        Quaternion_t Q_axis_z_A;

        Vector3d tmp_A1, tmp_A2;
        tmp_A1 << x2, 0, 0;
        tmp_A2 << 0, 0, 0;
        axis_z << 0, 0, 1;
        double angle_tmp_back = angle_tmp;
        NormalizeRotation(angle_tmp_back, axis_z, Q_axis_z_A);
        Rotation(tmp_A1, Q_axis_z_A, tmp_A2);

        tmp_A2(0) += x0;
        tmp_A2(1) += y0;
        tmp_A2(2) = 0;
        C = tmp_A2;
    }
    else if (abs(x0 - x1) < 0.0001 && abs(y0 - y1) > 0.0001) // x values are the same
    {
        double L = y1 - y0;
        double delta_y = (xi - X1) / (X2 - X1) * L;
        y2 = y0 + (0 + 1) * delta_y;
        x2 = x0;
        Vector3d tmp_u;
        tmp_u << x2, y2, 0;
        C = tmp_u;
    }
    else
    {
        cout << "Cannot transform (xi) to (xi, eta)!!\n";
        exit(0);
    }
    //cout << C.transpose() << endl;
};

bool comp(std::pair<size_t, double> &a, std::pair<size_t, double> &b)
{
    return a.second < b.second;
}

void Sort_trace_intersection_pnts(std::vector<Vector3d> &Intersec_traces_sub,
                                  const Vector3d A)
{
    if (Intersec_traces_sub.size() != 0 && Intersec_traces_sub.size() != 1)
    {
        std::vector<double> tmp_dis_1(Intersec_traces_sub.size());
        for (size_t k = 0; k < Intersec_traces_sub.size(); ++k)
        {
            tmp_dis_1[k] = pow((pow(Intersec_traces_sub[k](0) - A(0), 2) + pow(Intersec_traces_sub[k](1) - A(1), 2)), 0.5);
        }
        std::vector<pair<size_t, double>> tmp_dis_2(tmp_dis_1.size());
        for (size_t k = 0; k < tmp_dis_2.size(); ++k)
        {
            tmp_dis_2[k] = std::make_pair(k, tmp_dis_1[k]);
        }
        sort(tmp_dis_2.begin(), tmp_dis_2.end(), comp);

        std::vector<Vector3d> tmp_tt(Intersec_traces_sub.size());
        int tmp_u = 0;
        for (size_t k = 0; k < Intersec_traces_sub.size(); ++k)
        {
            if (k != 0)
            {
                if (abs(Intersec_traces_sub[tmp_dis_2[k].first](0) - tmp_tt[k - 1](0)) < 0.0001 && abs(Intersec_traces_sub[tmp_dis_2[k].first](1) - tmp_tt[k - 1](1)) < 0.0001)
                {
                    tmp_u++;
                }
                else
                {
                    tmp_tt[k - tmp_u] = Intersec_traces_sub[tmp_dis_2[k].first];
                    /*std::cout << "Intersec_traces[" << j << "][" << *it << "]"
                                  << ": " << Intersec_traces_sub[*it] << ";\n";*/
                }
            }
            else if (k == 0)
            {
                if (abs(Intersec_traces_sub[tmp_dis_2[k].first](0) - A(0)) < 0.0001 && abs(Intersec_traces_sub[tmp_dis_2[k].first](1) - A(1)) < 0.0001)
                {
                    tmp_u++;
                }
                else
                {
                    tmp_tt[k] = Intersec_traces_sub[tmp_dis_2[k].first];
                    /*std::cout << "Intersec_traces[" << j << "][" << *it << "]"
                                  << ": " << Intersec_traces_sub[*it] << ";\n";*/
                }
            }
        }
        Intersec_traces_sub.resize(0);
        Intersec_traces_sub.resize(tmp_tt.size() - tmp_u);
        for (size_t k = 0; k < Intersec_traces_sub.size(); ++k)
        {
            Intersec_traces_sub[k] = tmp_tt[k];
        }
        //std::cout << "found many intersection traces: \n";
    }
};

double Length_2d_segment(const Vector3d A, const Vector3d B)
{
    Vector3d C = A - B;
    return pow(pow(C(0), 2.) + pow(C(1), 2.), .5);
};

double unkown_vector_operator(const Vector3d p1, const Vector3d p2)
{
    return (p1(0) * p2(1) - p2(0) * p1(1));
}

bool pointInConvexPolygon(const Vector3d p, const std::vector<Vector3d> polygon)
{
    int i, iNext, i2Next;
    double preCross, nextCross;
    Vector3d v1, v2, v3;
    int polySize = polygon.size();

    if (polySize < 3)
    {
        return false;
    }
    for (i = 0; i < polySize; i++)
    {
        iNext = (i + 1) % polySize;
        i2Next = (iNext + 1) % polySize;

        v1 << polygon[i](0) - p(0), polygon[i](1) - p(1), 0;
        v2 << polygon[iNext](0) - p(0), polygon[iNext](1) - p(1), 0;
        preCross = unkown_vector_operator(v1, v2);

        v3 << polygon[i2Next](0) - p(0), polygon[i2Next](1) - p(1), 0;
        nextCross = unkown_vector_operator(v2, v3);

        if (preCross * nextCross < 0)
        {
            return false;
        }
    }

    return true;
}

bool Determine_IF_a_pnt_lies_in_a_horizontal_or_vertical_3D_square(const Vector3d A, const Vector3d B, const Vector3d C, const Vector3d D, const Vector3d PNT)
{
    //the center of this square
    Vector3d Center;
    Center = 0.5 * (A + C);

    Vector3d _A, _B, _C, _D, _PNT, _Normal;
    _A = A - Center;
    _B = B - Center;
    _C = C - Center;
    _D = D - Center;
    _PNT = PNT;
    _Normal = _A.cross(_B);
    //cout << "Normal: " << _Normal << endl;
    //cout << "Center: " << Center << endl;
    if (abs(_Normal(0)) < 0.0001 && abs(_Normal(1)) < 0.0001 && abs(_Normal(2)) > 0.0001) // it's a Top-zmax or bottom-zmin face
    {
        if (abs(_PNT(2) - Center(2)) < 0.001)
            return true;
        else
            return false;
    }
    else if (abs(_Normal(0)) < 0.0001 && abs(_Normal(1)) > 0.0001 && abs(_Normal(2)) < 0.0001) // it's a front-ymin or back-ymax face
    {
        if (abs(_PNT(1) - Center(1)) < 0.001)
            return true;
        else
            return false;
    }
    else if (abs(_Normal(0)) > 0.0001 && abs(_Normal(1)) < 0.0001 && abs(_Normal(2)) < 0.0001) // it's a left-xmin or right-xmax face
    {
        if (abs(_PNT(0) - Center(0)) < 0.001)
            return true;
        else
            return false;
    }
    else
    {
        cout << "Error! Model face is unusual!\n";
        exit(0);
    }
};

bool Determine_IF_a_pnt_lies_in_a_2D_line_segment(const Vector3d A, const Vector3d B, const Vector3d C)
{
    Vector3d _B, _C;
    _B = B - A;
    _C = C - A;

    double x_b = _B(0), y_b = _B(1);
    Vector3d axis_z;
    axis_z << 0, 0, -1;
    double R_angle_temp1 = atan2(y_b, x_b);

    Quaternion_t Q_axis_1;
    NormalizeRotation(R_angle_temp1, axis_z, Q_axis_1);

    Vector3d temp1;
    temp1 << x_b, y_b, 0;
    Vector3d temp2;
    Rotation(temp1, Q_axis_1, temp2);
    _B(0) = temp2(0);
    _B(1) = temp2(1);

    Vector3d temp3;
    temp3 = _C;
    Vector3d temp4;
    Rotation(temp3, Q_axis_1, temp4);
    _C(0) = temp4(0);
    _C(1) = temp4(1);

    if (abs(_C(1)) > 0.001) // y-axis
        return false;
    else
    {
        if (_C(0) > -0.001 && _C(0) - _B(0) < 0.001)
            return true;
        else
            return false;
    }
};

double Numerical_integration_linear_3(VectorXd Vec_q_normal, double L)
{
    if (Vec_q_normal.size() != 3)
    {
        throw Error_throw_ignore("Error! The boundary edge of an element should be divided into two sub-segments\n");
    }
    double q0 = Vec_q_normal(0);
    double q1 = Vec_q_normal(1);
    double q2 = Vec_q_normal(2);
    return (L * 0.5 * (q0 + 4 * q1 + q2)) / 3;
};

double Numerical_integration_linear_2(VectorXd Vec_q_normal, double L)
{
    if (Vec_q_normal.size() != 3)
    {
        throw Error_throw_ignore("Error! The boundary edge of an element should be divided into two sub-segments\n");
    }
    double q0 = Vec_q_normal(0);
    //double q1 = Vec_q_normal(1);
    double q2 = Vec_q_normal(2);

    return (L * (q0 + q2)) / 2;
};

double Heron_formula(Vector3d pnt0, Vector3d pnt2, Vector3d pnt4)
{
    double edge1, edge2, edge3;
    edge1 = (pnt0 - pnt2).norm();
    edge2 = (pnt2 - pnt4).norm();
    edge3 = (pnt0 - pnt4).norm();
    double p = 0.5 * (edge1 + edge2 + edge3);
    double A_area_2 = (p * (p - edge1) * (p - edge2) * (p - edge3));
    A_area_2 = pow(A_area_2, 0.5);
    return A_area_2;
};

void intersection_between_two_overlaped_line_segments(std::vector<Vector3d> Intersection_infinite,
                                                      Vector3d A,
                                                      Vector3d B,
                                                      std::vector<Vector3d> &Intersection_2)
{
    if (Intersection_2.size() != 2)
    {
        cout << "Found overlapped lines, the intersection vector should be re-initialized!\n";
        exit(0);
    }

    double x0, y0, x1, y1, x2, y2, x3, y3;
    x0 = Intersection_infinite[0](0);
    y0 = Intersection_infinite[0](1);
    x1 = Intersection_infinite[1](0);
    y1 = Intersection_infinite[1](1);
    x2 = A(0);
    y2 = A(1);
    x3 = B(0);
    y3 = B(1);

    //move to origin
    x1 -= x0;
    y1 -= y0;
    x2 -= x0;
    y2 -= y0;
    x3 -= x0;
    y3 -= y0;

    Vector3d axis_z;
    axis_z << 0, 0, -1;
    double R_angle_temp1 = atan2(y1, x1);

    Quaternion_t Q_axis_1;
    NormalizeRotation(R_angle_temp1, axis_z, Q_axis_1);

    Vector3d temp1;
    temp1 << x1, y1, 0;
    Vector3d temp2;
    Rotation(temp1, Q_axis_1, temp2);
    x1 = round(temp2(0), 1);
    y1 = round(temp2(1), 1);

    Vector3d temp3, temp4;
    temp3 << x2, y2, 0;
    Rotation(temp3, Q_axis_1, temp4);
    x2 = round(temp4(0), 1);
    y2 = round(temp4(1), 1);

    Vector3d temp5, temp6;
    temp5 << x3, y3, 0;
    Rotation(temp5, Q_axis_1, temp6);
    x3 = round(temp6(0), 1);
    y3 = round(temp6(1), 1);

    if (abs(y1) > 0.001 || abs(y2) > 0.001 || abs(y3) > 0.001)
    {
        cout << y0 << ", " << y1 << ", " << y2 << ", " << y3 << "\n";
        cout << "line 1: \n";
        cout << Intersection_infinite[0].transpose() << endl;
        cout << Intersection_infinite[1].transpose() << endl;
        cout << "line 2: \n";
        cout << A.transpose() << endl;
        cout << B.transpose() << endl;
        cout << "rotation fails!\n";
        exit(0);
    }
    Vector3d A1;
    Vector3d B1;
    Vector3d C1;
    Vector3d D1;
    A1 << 0, 0, 0;
    B1 << x1, y1, 0;
    C1 << x2, y2, 0;
    D1 << x3, y3, 0;

    std::vector<Vector3d> Intersection_4 = {A1, B1};
    std::vector<Vector3d> Intersection_5 = {C1, D1};
    Intersection_of_1D_intervals(Intersection_4, Intersection_5, Intersection_2);

    //rotation back and
    //move back
    NormalizeRotation(-R_angle_temp1, axis_z, Q_axis_1);

    //cout << "ppp\n";
    Vector3d CRTR2;
    CRTR2 << 0, 0, 0;
    Rotation(Intersection_2[0], Q_axis_1, CRTR2);
    Intersection_2[0] << round(CRTR2(0), 4) + x0, round(CRTR2(1), 4) + y0, 0;

    CRTR2 << 0, 0, 0;
    Rotation(Intersection_2[1], Q_axis_1, CRTR2);
    Intersection_2[1] << round(CRTR2(0), 4) + x0, round(CRTR2(1), 4) + y0, 0;
    //cout << "pppoo\n";
};

string To_string_with_width(size_t val, size_t width)
{
    std::ostringstream oss;
    oss.width(width);
    oss.fill('0');
    oss << val;
    return oss.str();
}

void time_counter_start(Time_Pnt &start)
{
    start = std::chrono::steady_clock::now();
};

void time_counter_end(Time_Pnt start,
                      Time_Pnt &end,
                      string content,
                      string unit)
{
    end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::micro> elapsed = end - start; // std::micro time (us)

    if (unit == "in_hours")
        std::cout << "Running time of " << content << ": " << (((double)(elapsed.count() * 1.0) * (0.000001)) / 60.00) / 60.00 << " hours" << std::endl;
    else if (unit == "in_minutes")
        std::cout << "Running time of " << content << ": " << (((double)(elapsed.count() * 1.0) * (0.000001)) / 60.00) << " minutes" << std::endl;
    else if (unit == "in_seconds")
        std::cout << "Running time of " << content << ": " << (((double)(elapsed.count() * 1.0) * (0.000001))) << " seconds" << std::endl;
    else
    {
        string AS = "Please define the unit with input 'in_hours', 'in_minutes', or 'in_seconds'!\n";
        throw Error_throw_pause(AS);
    }
}