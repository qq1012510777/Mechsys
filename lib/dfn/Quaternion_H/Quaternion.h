#pragma once
#include "Dense"

// Std Lib
#include <cmath>
#include <iostream>
using namespace std;

using namespace Eigen;

typedef Eigen::Vector4d Quaternion_t;

inline void NormalizeRotation(double Theta, Vector3d const &Axis, Quaternion_t &C)
{
    if (Axis.norm() < 1.0e-12)
    {
        std::cout << "Quaternion: the norm of the axis is too small, please chose a different one\n";
        exit(0);
    }
    Vector3d A = Axis / Axis.norm();

    C(0) = cos(Theta / 2.0);
    C(1) = A(0) * sin(Theta / 2.0);
    C(2) = A(1) * sin(Theta / 2.0);
    C(3) = A(2) * sin(Theta / 2.0);
}

inline void Conjugate_Q(Quaternion_t const &A, Quaternion_t &C)
{
    C(0) = A(0);
    C(1) = -A(1);
    C(2) = -A(2);
    C(3) = -A(3);
}

inline void GetVector(Quaternion_t const &A, Vector3d &C)
{
    C(0) = A(1);
    C(1) = A(2);
    C(2) = A(3);
}

inline void SetQuaternion(double Scalar, Vector3d const &A, Quaternion_t &C)
{
    C(0) = Scalar;
    C(1) = A(0);
    C(2) = A(1);
    C(3) = A(2);
}

inline void QuaternionProduct(Quaternion_t const &A, Quaternion_t const &B, Quaternion_t &C)
{
    Vector3d t1, t2;
    GetVector(A, t1);
    GetVector(B, t2);
    double scalar = A(0) * B(0) - t1.dot(t2);
    Vector3d vector = A(0) * t2 + B(0) * t1 + t1.cross(t2);
    SetQuaternion(scalar, vector, C);
}

inline void Rotation(Vector3d const &A, Quaternion_t const &B, Vector3d &C)
{
    Quaternion_t t1, t2, t3;
    SetQuaternion(0.0, A, t1);
    QuaternionProduct(B, t1, t2);
    Conjugate_Q(B, t3);
    QuaternionProduct(t2, t3, t1);
    GetVector(t1, C);
}
