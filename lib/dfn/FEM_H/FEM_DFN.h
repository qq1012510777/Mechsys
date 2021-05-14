#pragma once

#include "../Mesh_H/Mesh_DFN.h"

namespace DFN
{

class FEM_DFN
{
public:
    FEM_DFN(DFN::Mesh_DFN DFN_mesh);
    void Assemble_overall_matrix(DFN::Mesh_DFN DFN_mesh);
};

inline FEM_DFN::FEM_DFN(DFN::Mesh_DFN DFN_mesh){

};

inline void FEM_DFN::Assemble_overall_matrix(DFN::Mesh_DFN DFN_mesh)
{
    for (size_t i = 0; i < DFN_mesh.JM.size(); ++i)
    {
        for (size_t j = 0; j < DFN_mesh.JM[i].size(); ++j)
        {
            Eigen::RowVectorXd w, xi, eta;
            w = Eigen::RowVectorXd::Zero(6);
            w << 0.1713244923791700, 0.3607615730481380, 0.4679139345726910, 0.1713244923791700, 0.3607615730481380, 0.4679139345726910;
            xi = Eigen::RowVectorXd::Zero(6);
            xi << 0.9324695142031520, 0.6612093864662640, 0.2386191860831960, -0.9324695142031520, -0.6612093864662640, -0.2386191860831960;
            eta = xi;

            MatrixXd D1e;
            D1e = MatrixXd::Zero(6, 6);

            for (size_t ik = 0; ik < (size_t)D1e.rows(); ++ik)
            {
                for (size_t jk = 0; jk < (size_t)D1e.cols(); ++jk)
                {
                    Eigen::MatrixXd JXYe_x, JXYe_y;
                    JXYe_x = Eigen::MatrixXd::Zero(6, 1);
                    JXYe_y = Eigen::MatrixXd::Zero(6, 1);

                    for (size_t kk = 0; kk < (size_t)JXYe_x.rows(); ++kk)
                    {
                        JXYe_x(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](kk)](0);
                        JXYe_y(kk, 0) = DFN_mesh.JXY[i][DFN_mesh.JM[i][j](kk)](1);
                    }

                    VectorXd pd_N_over_pd_x, pd_N_over_pd_y;
                    pd_N_over_pd_x = VectorXd::Zero(6);
                    pd_N_over_pd_y = VectorXd::Zero(6);

                    MatrixXd Jacobi;
                    Jacobi = MatrixXd::Zero(2, 2);

                    Eigen::VectorXd PHI;
                    PHI = Eigen::VectorXd::Zero(6);

                    double Phi_1;
                    double Phi_2;
                    double Phi_3;
                    double Phi_4;
                    double Phi_5;
                    double Phi_6;
                    
                    this->PHI_shape_function(xi(ik),
                                       eta(jk),
                                       Phi_1,
                                       Phi_2,
                                       Phi_3,
                                       Phi_4,
                                       Phi_5,
                                       Phi_6);
                    PHI << Phi_1, Phi_4, Phi_2, Phi_5, Phi_3, Phi_6;

                    p_PHI_over_p_x_and_y(xi(ik), eta(jk), JXYe_x, JXYe_y, pd_N_over_pd_x, pd_N_over_pd_y, Jacobi);

                    D1e = D1e + 1. / Kper * w(ik) * w(jk) * Jacobi.determinant() * PHI * PHI.transpose();
                };
            };
        }
    }
};

}; // namespace DFN
