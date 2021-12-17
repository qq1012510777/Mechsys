#pragma once
#include "../Mesh_DFN_linear.h"

namespace DFN
{
class SEP_edge_numbering
{
public:
    vector<MatrixXs> Sep_edges_NO_of_ele_frac;
    size_t Num_edges_sep = 0;

public:
    SEP_edge_numbering(const DFN::Mesh_DFN_linear mesh);
    void Matlab_plot(string FileKey_mat,
                     string FileKey_m,
                     DFN::Mesh_DFN_linear mesh);
};

inline SEP_edge_numbering::SEP_edge_numbering(const DFN::Mesh_DFN_linear mesh)
{
    Sep_edges_NO_of_ele_frac.resize(mesh.element_2D.size());

    size_t edgeNO = 1;
    for (size_t i = 0; i < Sep_edges_NO_of_ele_frac.size(); ++i)
    {
        Sep_edges_NO_of_ele_frac[i].resize(mesh.element_2D[i].rows(), mesh.element_2D[i].cols());

        for (size_t j = 0; j < (size_t)mesh.element_2D[i].rows(); ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                Sep_edges_NO_of_ele_frac[i](j, k) = edgeNO;
                edgeNO++;
            }
        }

    }
    Num_edges_sep = edgeNO - 1;
};

void SEP_edge_numbering::Matlab_plot(string FileKey_mat,
                                     string FileKey_m,
                                     DFN::Mesh_DFN_linear mesh)
{
    const char *filename = FileKey_mat.c_str();

    DFN::MATLAB_DATA_API M1_;
    M1_.Write_mat(filename, "w", 1, 1, 1, {0}, "nothing_");

    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i)
    {
        //cout << i << endl;
        MatrixXs ele_2D_frac = mesh.element_2D[i];

        vector<double> pData1(ele_2D_frac.rows() * 3);

        for (size_t j = 0; j < (size_t)ele_2D_frac.rows() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / ele_2D_frac.rows()); // column
            l = j % ele_2D_frac.rows();       // row

            pData1[j] = ele_2D_frac(l, k);
        }
        M1_.Write_mat(filename, "u", ele_2D_frac.rows() * 3,
                      ele_2D_frac.rows(), 3, pData1,
                      "element_2D_Frac_" + to_string(i + 1));

        //----------------------------------
        vector<double> pData2(mesh.coordinate_2D[i].rows() * 2);

        for (size_t j = 0; j < (size_t)mesh.coordinate_2D[i].rows() * 2; ++j)
        {
            size_t k, l;
            k = ceil(j / mesh.coordinate_2D[i].rows()); // column
            l = j % mesh.coordinate_2D[i].rows();       // row

            pData2[j] = mesh.coordinate_2D[i].coeffRef(l, k);
        }
        M1_.Write_mat(filename, "u", mesh.coordinate_2D[i].rows() * 2,
                      mesh.coordinate_2D[i].rows(), 2, pData2,
                      "coordinate_2D_Frac_" + to_string(i + 1));

        //----------------------------------
        vector<double> pData3(this->Sep_edges_NO_of_ele_frac[i].rows() * 3);

        for (size_t j = 0; j < (size_t)this->Sep_edges_NO_of_ele_frac[i].rows() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / this->Sep_edges_NO_of_ele_frac[i].rows()); // column
            l = j % this->Sep_edges_NO_of_ele_frac[i].rows();       // row

            pData3[j] = this->Sep_edges_NO_of_ele_frac[i](l, k);
            //cout << pData3[j] << endl;
        }
        M1_.Write_mat(filename, "u",
                      this->Sep_edges_NO_of_ele_frac[i].rows() * 3,
                      this->Sep_edges_NO_of_ele_frac[i].rows(),
                      3,
                      pData3,
                      "Sep_edges_NO_of_ele_frac_" + to_string(i + 1));
    }
};
}; // namespace DFN