#pragma once
#include "GLOB_edge_numbering.h"
#include "SEP_edge_numbering.h"

namespace DFN
{
class MAT_plot_SEP_GLOB_edges
{
public:
    MAT_plot_SEP_GLOB_edges(string FileKey_mat,
                            string FileKey_m,
                            DFN::Mesh_DFN_linear mesh,
                            DFN::SEP_edge_numbering sep,
                            DFN::GLOB_edge_numbering glob);
};

MAT_plot_SEP_GLOB_edges::MAT_plot_SEP_GLOB_edges(string FileKey_mat,
                                                 string FileKey_m,
                                                 DFN::Mesh_DFN_linear mesh,
                                                 DFN::SEP_edge_numbering sep,
                                                 DFN::GLOB_edge_numbering glob)
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
        vector<double> pData3(sep.Sep_edges_NO_of_ele_frac[i].rows() * 3);

        for (size_t j = 0; j < (size_t)sep.Sep_edges_NO_of_ele_frac[i].rows() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / sep.Sep_edges_NO_of_ele_frac[i].rows()); // column
            l = j % sep.Sep_edges_NO_of_ele_frac[i].rows();       // row

            pData3[j] = sep.Sep_edges_NO_of_ele_frac[i](l, k);
            //cout << pData3[j] << endl;
        }
        //cout << sep.Sep_edges_NO_of_ele_frac[0].row(0) << endl;
        M1_.Write_mat(filename, "u",
                      sep.Sep_edges_NO_of_ele_frac[i].rows() * 3,
                      sep.Sep_edges_NO_of_ele_frac[i].rows(),
                      3,
                      pData3,
                      "Sep_edges_NO_of_ele_frac_" + to_string(i + 1));
        //----------------------------------
        vector<double> pData4(glob.GLOB_edges_NO_of_ele_frac[i].rows() * 3);

        for (size_t j = 0; j < (size_t)glob.GLOB_edges_NO_of_ele_frac[i].rows() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / glob.GLOB_edges_NO_of_ele_frac[i].rows()); // column
            l = j % glob.GLOB_edges_NO_of_ele_frac[i].rows();       // row

            pData4[j] = glob.GLOB_edges_NO_of_ele_frac[i](l, k);
        }

        M1_.Write_mat(filename, "u",
                      glob.GLOB_edges_NO_of_ele_frac[i].rows() * 3,
                      glob.GLOB_edges_NO_of_ele_frac[i].rows(),
                      3,
                      pData4,
                      "GLOB_edges_NO_of_ele_frac_" + to_string(i + 1));
    }

    vector<double> pData4(mesh.coordinate_3D.rows() * 3);
    vector<double> pData5(mesh.element_3D.rows() * 3);

    for (size_t j = 0; j < (size_t)mesh.coordinate_3D.rows() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / mesh.coordinate_3D.rows()); // column
        l = j % mesh.coordinate_3D.rows();       // row

        pData4[j] = mesh.coordinate_3D(l, k);
        //cout << pData4[j] << endl;
    }

    for (size_t j = 0; j < (size_t)mesh.element_3D.rows() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / mesh.element_3D.rows()); // column
        l = j % mesh.element_3D.rows();       // row

        pData5[j] = mesh.element_3D(l, k);
    }

    string FracJXY3D_s = "coordinate_3D";
    string FracJM_s = "element_3D";

    M1_.Write_mat(filename, "u", mesh.coordinate_3D.rows() * 3,
                  mesh.coordinate_3D.rows(), 3, pData4, FracJXY3D_s);
    M1_.Write_mat(filename, "u", mesh.element_3D.rows() * 3,
                  mesh.element_3D.rows(), 3, pData5, FracJM_s);

    // m file
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;";
    oss << "load('" << FileKey_mat << "');\n";

    size_t figure1 = 1;

    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i)
    {
        oss << "figure(" << figure1 << "); title('SEP. edge NO FRAC " << figure1 << "')\n";
        oss << "patch('Vertices', coordinate_2D_Frac_" << figure1 << ", 'Faces', element_2D_Frac_" << figure1 << ", 'FaceVertexCData', zeros(size(coordinate_2D_Frac_1, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);\n";
        oss << "for i = 1:size(element_2D_Frac_" << figure1 << ", 1)\n";
        oss << "\tnode1 = element_2D_Frac_" << figure1 << "(i, [1, 2, 3]);\n";
        oss << "\tnode2 = element_2D_Frac_" << figure1 << "(i, [2, 3, 1]);\n";
        oss << "\tcoord = coordinate_2D_Frac_" << figure1 << "(node1, :) + coordinate_2D_Frac_" << figure1 << "(node2, :);\n";
        oss << "\tcoord = coord .* 0.5;\n";
        oss << "\ttext(coord(:, 1), coord(:, 2), num2str([Sep_edges_NO_of_ele_frac_" << figure1 << "(i, :)]')); hold on\n";

        oss << "\tfor j = 1:3\n";
        oss << "\t\tnode1 = element_2D_Frac_" << figure1 << "(i, j);\n";
        oss << "\t\ttext(coordinate_2D_Frac_" << figure1 << "(node1, 1), coordinate_2D_Frac_" << figure1 << "(node1, 2), num2str(node1), 'color', 'red'); hold on\n";
        oss << "\tend\n";
        oss << "end\n";
        oss << "\n";
        figure1++;
    }
    figure1--;
    size_t ysu = figure1;

    figure1 = 1;
    for (size_t i = 0; i < mesh.Frac_Tag.size(); ++i)
    {
        oss << "figure(" << ysu + 1 << "); title('GLOB. edge NO FRAC " << figure1 << "')\n";
        oss << "patch('Vertices', coordinate_3D, 'Faces', element_2D_Frac_" << figure1 << ", 'FaceVertexCData', zeros(size(coordinate_3D, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0); hold on; view(3)\n";
        oss << "for i = 1:size(element_2D_Frac_" << figure1 << ", 1)\n";
        oss << "\tnode1 = element_2D_Frac_" << figure1 << "(i, [1, 2, 3]);\n";
        oss << "\tnode2 = element_2D_Frac_" << figure1 << "(i, [2, 3, 1]);\n";
        oss << "\tcoord = coordinate_3D(node1, :) + coordinate_3D(node2, :);\n";
        oss << "\tcoord = coord .* 0.5;\n";
        oss << "\ttext(coord(:, 1), coord(:, 2), coord(:, 3), num2str([GLOB_edges_NO_of_ele_frac_" << figure1 << "(i, :)]')); hold on;\n";
        oss << "end\n";
        oss << "\n";
        figure1++;
    }

    oss.close();
}
}; // namespace DFN