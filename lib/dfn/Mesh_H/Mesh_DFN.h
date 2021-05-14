#pragma once
#include "../DFN_H/Domain_WL.h"
#include "../DFN_H/Remove_unnecess_frac.h"
#include "../Geometry_H/Handle_Polygon_2D_with_extream_traces.h"
#include "../Geometry_H/Line_seg_2D.h"
#include "../Geometry_H/Polygon_convex_2D.h"
#include "../Geometry_H/Polygon_convex_2D_with_traces.h"
#include "../Geometry_H/Polygon_convex_3D.h"
#include "../Geometry_H/Splitting_Polygon_convex_2D_with_traces.h"
#include "../Geometry_H/Vector_2.h"
#include "Mesh_Polygon_2D_with_traces.h"
#include "mat.h"

typedef struct If_Model_bound_pnt
{
    bool If_model_top = false;
    bool If_model_bottom = false;
    bool If_model_front = false;
    bool If_model_back = false;
    bool If_model_left = false;
    bool If_model_right = false;
    //Top-zmax, bottom-zmin, front-ymin, back-ymax, left-xmin, right-xmax
} Mo_B_pnt;

namespace DFN
{
class Mesh_DFN
{
public:
    std::vector<size_t> Frac_Tag;
    std::vector<std::vector<spe_pnt>> Pnt_attribute;   // frac bond, trace
    std::vector<std::vector<Mo_B_pnt>> Pnt_attri_Mo_B; // Model bound

    std::vector<std::vector<Vector2d>> JXY;
    ///< JXY.size() is the number of fractures,
    ///< JXY[0].size() is the number of vertices of

    std::vector<std::vector<Vector3d>> JXY_3D;

    std::vector<std::vector<RowVector6i>> JM;
    ///< each Vector3d is the ID values of
    ///< the three nodes

    std::vector<std::vector<size_t>> Trace_Tag; // record point Tags belonging to traces of each fracture
    std::vector<std::vector<Eigen::Vector3d>> Guide_frac;
    ///< (0): an index referring if this is
    //a repetitive point;
    //(1): i; (2): j;(if repetitive)

    size_t FEM_matrix_dimension;

public:
    Mesh_DFN(DFN::Domain dom);
    void Create_structure_of_a_frac(DFN::Domain dom, const size_t FracID, std::vector<Vector3d> &Verts, vector<std::pair<Vector3d, Vector3d>> &Traces, DFN::Remove_unnecess_frac Remove_f);
    void Matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom);
    void Ident_Model_Bound_Pnt(DFN::Domain dom);
    void Find_repetitive_pnt(DFN::Domain dom);
};

inline Mesh_DFN::Mesh_DFN(DFN::Domain dom)
{
    DFN::Remove_unnecess_frac Remove_f(dom);

    for (size_t i = 0; i < Remove_f.Listofclusters_remove_fracs.size(); ++i)
    {
        for (size_t j = 0; j < Remove_f.Listofclusters_remove_fracs[i].size(); ++j)
        {
            if (Remove_f.Listofclusters_remove_fracs[i][j] != -1)
            {

                this->Frac_Tag.push_back(Remove_f.Listofclusters_remove_fracs[i][j]);
                size_t FracID = Remove_f.Listofclusters_remove_fracs[i][j];
                std::vector<Vector3d> Verts_frac;
                std::vector<std::pair<Vector3d, Vector3d>> Traces_frac;
                this->Create_structure_of_a_frac(dom, FracID, Verts_frac, Traces_frac, Remove_f);

                DFN::Polygon_convex_3D poly{Verts_frac};
                Vector3d temp1;
                DFN::Vector_2 v(poly.Normal_vector, temp1);

                double R_angle_temp1 = 0;
                double x_temp = poly.Beta;
                R_angle_temp1 = x_temp * M_PI / 180;
                Quaternion_t Q_axis_1;

                std::vector<Vector3d> Verts_frac_rotated, Traces_frac_f, Traces_frac_f_rotated;
                Traces_frac_f.resize(Traces_frac.size() * 2);
                for (size_t k = 0; k < Traces_frac_f.size(); k = k + 2)
                {
                    size_t pd = k / 2;
                    Traces_frac_f[k] = Traces_frac[pd].first;
                    Traces_frac_f[k + 1] = Traces_frac[pd].second;
                }

                if (poly.Beta > 0.0001)
                {
                    //cout << "1;\n";
                    Verts_frac_rotated.resize(Verts_frac.size());
                    Traces_frac_f_rotated.resize(Traces_frac_f.size());
                    DFN::Rotation_verts R2(Verts_frac, R_angle_temp1, Q_axis_1, temp1, Verts_frac_rotated);
                    DFN::Rotation_verts R3(Traces_frac_f, R_angle_temp1, Q_axis_1, temp1, Traces_frac_f_rotated);
                }
                else
                {
                    Verts_frac_rotated = Verts_frac;
                    Traces_frac_f_rotated = Traces_frac_f;
                }

                //move to x-y plane
                double z_coord = Verts_frac_rotated[0](2);

                for (size_t k = 0; k < Verts_frac_rotated.size(); ++k)
                    Verts_frac_rotated[k](2) = 0;

                for (size_t k = 0; k < Traces_frac_f_rotated.size(); ++k)
                    Traces_frac_f_rotated[k](2) = 0;

                std::vector<DFN::Line_seg_2D> Traces_Line(Traces_frac.size());

                for (size_t k = 0; k < Traces_frac_f_rotated.size(); k = k + 2)
                {
                    size_t pd = k / 2;
                    Traces_Line[pd].Re_constructor(Vector2d{Traces_frac_f_rotated[k](0), Traces_frac_f_rotated[k](1)}, Vector2d{Traces_frac_f_rotated[k + 1](0), Traces_frac_f_rotated[k + 1](1)});
                }
                // there should be a pre_handle code to divide polygon2D if extream traces exist
                // extream means: the trace nearly intersect edges but actually not
                // and very tiny trace but not a point

                DFN::Polygon_convex_2D polyframe{Verts_frac_rotated};

                DFN::Polygon_convex_2D_with_traces poly2D{polyframe, Traces_Line};

                DFN::Splitting_Polygon_convex_2D_with_traces splitting_poly2D{poly2D, 1};

                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> neigh_shared_A;

                size_t NO_Nodes_p_A = 0;

                DFN::Mesh_Polygon_2D_with_traces mesh{
                    splitting_poly2D.Pnt_sets,
                    splitting_poly2D.Topo_struct,
                    29,
                    0.001,
                    3,
                    neigh_shared_A,
                    NO_Nodes_p_A};

                mesh.Identifying_Frac_bound(polyframe, Traces_Line);
                this->Pnt_attribute.push_back(mesh.Pnt_attribute);

                this->JXY.push_back(mesh.JXY);
                this->JM.push_back(mesh.JM);
                this->Trace_Tag.push_back(mesh.Trace_Tag);

                std::vector<Vector3d> Verts_2D_tmp(mesh.JXY.size());
                for (size_t k = 0; k < Verts_2D_tmp.size(); ++k)
                    Verts_2D_tmp[k] << mesh.JXY[k](0), mesh.JXY[k](1), z_coord; // move back
                std::vector<Vector3d> Verts_rotated_back(mesh.JXY.size());

                if (poly.Beta > 0.0001)
                {
                    DFN::Rotation_verts R2(Verts_2D_tmp, -R_angle_temp1, Q_axis_1, temp1, Verts_rotated_back);
                }
                else
                {
                    Verts_rotated_back = Verts_2D_tmp;
                }
                this->JXY_3D.push_back(Verts_rotated_back);
            }
        }
    }

    this->Ident_Model_Bound_Pnt(dom);
    this->Find_repetitive_pnt(dom);
};

inline void Mesh_DFN::Create_structure_of_a_frac(DFN::Domain dom, const size_t FracID, std::vector<Vector3d> &Verts, std::vector<std::pair<Vector3d, Vector3d>> &Traces, DFN::Remove_unnecess_frac Remove_f)
{

    Verts.resize(dom.Fractures[FracID].Verts_trim.size());
    for (size_t i = 0; i < Verts.size(); ++i)
        Verts[i] = dom.Fractures[FracID].Verts_trim[i];

    std::set<size_t>::iterator it_a = dom.Fractures[FracID].Intersect_other_frac_after_trim.begin();

    while (it_a != dom.Fractures[FracID].Intersect_other_frac_after_trim.end())
    {
        size_t other_fracID = *it_a;

        // if other frac is removed?
        bool removed = false;
        for (size_t i = 0; i < Remove_f.Listofclusters_remove_fracs.size(); ++i)
        {
            std::vector<int>::iterator it_r = find(Remove_f.Listofclusters_remove_fracs[i].begin(), Remove_f.Listofclusters_remove_fracs[i].end(), other_fracID);
            if (it_r == Remove_f.Listofclusters_remove_fracs[i].end())
            {
                removed = true;
                break;
            }
        }
        if (removed == true)
        {
            it_a++;
            continue;
        }

        size_t ID_1 = FracID < other_fracID ? FracID : other_fracID;
        size_t ID_2 = FracID > other_fracID ? FracID : other_fracID;

        std::pair<size_t, size_t> Trace_key = std::make_pair(ID_1, ID_2);

        Traces.push_back(std::make_pair(dom.Intersections[Trace_key].first, dom.Intersections[Trace_key].second));
        it_a++;
    }
};

void Mesh_DFN::Matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom)
{
    // mat file
    //cout << this->JM[0][0].transpose() << endl;
    const char *filename = FileKey_mat.c_str();
    MATFile *pMatFile;
    pMatFile = matOpen(filename, "w");

    if (!pMatFile)
    {
        cout << "cannot create mat file in class Mesh_DFN\n";
        exit(0);
    }
    // frac
    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        //cout << "i: " << i << endl;
        size_t len = dom.Fractures[Frac_Tag[i]].Verts_trim.size(); // number of verts

        double *pData1;
        double *pData2;
        double *pData3;
        double *pData4;
        double *pData5;
        double *pData6;
        double *pData7;
        pData1 = (double *)mxCalloc(len, sizeof(double));
        pData2 = (double *)mxCalloc(len, sizeof(double));
        pData3 = (double *)mxCalloc(len, sizeof(double));
        pData4 = (double *)mxCalloc(this->JXY_3D[i].size() * 3, sizeof(double));
        pData5 = (double *)mxCalloc(this->JM[i].size() * 6, sizeof(double));
        pData6 = (double *)mxCalloc(this->JXY_3D[i].size() * 2, sizeof(double));
        pData7 = (double *)mxCalloc(this->JXY_3D[i].size() * 6, sizeof(double));

        mxArray *pMxArray1;
        mxArray *pMxArray2;
        mxArray *pMxArray3;
        mxArray *pMxArray4;
        mxArray *pMxArray5;
        mxArray *pMxArray6;
        mxArray *pMxArray7;

        pMxArray1 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray2 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray3 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray4 = mxCreateDoubleMatrix(this->JXY_3D[i].size(), 3, mxREAL);
        pMxArray5 = mxCreateDoubleMatrix(this->JM[i].size(), 6, mxREAL);
        pMxArray6 = mxCreateDoubleMatrix(this->JXY_3D[i].size(), 2, mxREAL);
        pMxArray7 = mxCreateDoubleMatrix(this->JXY_3D[i].size(), 6, mxREAL);

        if (!pMxArray1 || !pMxArray2 || !pMxArray3 || !pMxArray4 || !pMxArray5 || !pMxArray6 || !pMxArray7)
        {
            cout << "cannot create pMxArray in class Mesh_DFN\n";
            exit(0);
        }

        if (!pData1 || !pData2 || !pData3 || !pData4 || !pData5 || !pData6 || !pData7)
        {
            cout << "cannot create pData in class Mesh_DFN\n";
            exit(0);
        }

        for (size_t j = 0; j < len; j++)
        {
            pData1[j] = dom.Fractures[Frac_Tag[i]].Verts_trim[j](0);
            pData2[j] = dom.Fractures[Frac_Tag[i]].Verts_trim[j](1);
            pData3[j] = dom.Fractures[Frac_Tag[i]].Verts_trim[j](2);
        }

        for (size_t j = 0; j < this->JXY_3D[i].size() * 3; ++j)
        {
            size_t k, l;
            k = ceil(j / this->JXY_3D[i].size()); // column
            l = j % this->JXY_3D[i].size();       // row

            pData4[j] = this->JXY_3D[i][l](k);
        }

        for (size_t j = 0; j < this->JM[i].size() * 6; ++j)
        {
            size_t k, l;
            k = ceil(j / this->JM[i].size()); // column
            l = j % this->JM[i].size();       // row

            pData5[j] = this->JM[i][l](k) + 1;
        }

        for (size_t j = 0; j < this->JXY_3D[i].size() * 2; ++j)
        {
            size_t k, l;
            k = ceil(j / this->JXY_3D[i].size()); // column
            l = j % this->JXY_3D[i].size();       // row

            pData6[j] = 0;

            if (k == 0)
                if (this->Pnt_attribute[i][l].If_frac_bound == true)
                    pData6[j] = 1;

            if (k == 1)
                if (this->Pnt_attribute[i][l].If_trace == true)
                    pData6[j] = 1;
        }

        for (size_t j = 0; j < this->JXY_3D[i].size() * 6; ++j)
        {
            size_t k, l;
            k = ceil(j / this->JXY_3D[i].size()); // column
            l = j % this->JXY_3D[i].size();       // row

            pData7[j] = 0;

            if (k == 0)
                if (this->Pnt_attri_Mo_B[i][l].If_model_top == true)
                    pData7[j] = 1;
            if (k == 1)
                if (this->Pnt_attri_Mo_B[i][l].If_model_bottom == true)
                    pData7[j] = 1;
            if (k == 2)
                if (this->Pnt_attri_Mo_B[i][l].If_model_front == true)
                    pData7[j] = 1;
            if (k == 3)
                if (this->Pnt_attri_Mo_B[i][l].If_model_back == true)
                    pData7[j] = 1;
            if (k == 4)
                if (this->Pnt_attri_Mo_B[i][l].If_model_left == true)
                    pData7[j] = 1;
            if (k == 5)
                if (this->Pnt_attri_Mo_B[i][l].If_model_right == true)
                    pData7[j] = 1;
        }

        mxSetData(pMxArray1, pData1);
        mxSetData(pMxArray2, pData2);
        mxSetData(pMxArray3, pData3);
        mxSetData(pMxArray4, pData4);
        mxSetData(pMxArray5, pData5);
        mxSetData(pMxArray6, pData6);
        mxSetData(pMxArray7, pData7);

        string ft = to_string(i);

        string Fracx = "Frac_" + ft + "_x";
        string Fracy = "Frac_" + ft + "_y";
        string Fracz = "Frac_" + ft + "_z";
        string FracJXY3D = "Frac_" + ft + "_JXY3D";
        string FracJM = "Frac_" + ft + "_JM";
        string FracPnt_attribute = "Frac_" + ft + "_Pnt_attribute";
        string FracPnt_attri_Mo_B = "Frac_" + ft + "_Pnt_attri_Mo_B";

        const char *Fracx_s = Fracx.c_str();
        const char *Fracy_s = Fracy.c_str();
        const char *Fracz_s = Fracz.c_str();
        const char *FracJXY3D_s = FracJXY3D.c_str();
        const char *FracJM_s = FracJM.c_str();
        const char *FracPnt_attribute_s = FracPnt_attribute.c_str();
        const char *FracPnt_attri_Mo_B_s = FracPnt_attri_Mo_B.c_str();

        matPutVariable(pMatFile, Fracx_s, pMxArray1);
        matPutVariable(pMatFile, Fracy_s, pMxArray2);
        matPutVariable(pMatFile, Fracz_s, pMxArray3);
        matPutVariable(pMatFile, FracJXY3D_s, pMxArray4);
        matPutVariable(pMatFile, FracJM_s, pMxArray5);
        matPutVariable(pMatFile, FracPnt_attribute_s, pMxArray6);
        matPutVariable(pMatFile, FracPnt_attri_Mo_B_s, pMxArray7);

        mxFree(pData1);
        mxFree(pData2);
        mxFree(pData3);
        mxFree(pData4);
        mxFree(pData5);
        mxFree(pData6);
        mxFree(pData7);
    }

    matClose(pMatFile);

    // m file
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "load('" << FileKey_mat << "');\n";
    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        oss << "fill3([Frac_" << i << "_x; Frac_" << i << "_x(1,1)], [Frac_" << i << "_y; Frac_" << i << "_y(1,1)], [Frac_" << i << "_z; Frac_" << i << "_z(1,1)], [rand rand rand]);\nhold on;\n";
    }
    oss << "\n\n";
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << dom.Surfaces[i].Verts[j](0) << " " << dom.Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](1) << " " << dom.Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << dom.Surfaces[i].Verts[j](2) << " " << dom.Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on; hold on;\n";
        }
    }
    double xmin_1 = dom.Model_domain(4), xmax_1 = dom.Model_domain(5);
    double ymin_1 = dom.Model_domain(2), ymax_1 = dom.Model_domain(3);
    double zmin_1 = dom.Model_domain(1), zmax_1 = dom.Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

    oss << "\n\n";
    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        oss << "[m" << i << ", ~] = size(Frac_" << i << "_JXY3D);\n";
        oss << "P" << i << " = patch('Vertices', Frac_" << i << "_JXY3D, 'Faces', Frac_" << i << "_JM, 'FaceVertexCData', zeros(m" << i << ", 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);\n";
        oss << "hold on;\n";
    }

    oss << "\n\n\n";

    oss << "Frac_bond = [];\n";
    oss << "Frac_count = 1;\n";

    oss << "Trace_bond = [];\n";
    oss << "Trace_count = 1;\n";

    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        oss << "[m" << i << ", ~] = size(Frac_" << i << "_Pnt_attribute);\n";
        oss << "for i = 1 : m" << i << "\n";
        oss << "\tfor j = 1 : 2\n";
        oss << "\t\tif (Frac_" << i << "_Pnt_attribute(i,j) == 1 && j == 1)\n";
        oss << "\t\t\tFrac_bond(Frac_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tFrac_count = Frac_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n\n";

        oss << "\t\tif (Frac_" << i << "_Pnt_attribute(i,j) == 1 && j == 2)\n";
        oss << "\t\t\tTrace_bond(Trace_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tTrace_count = Trace_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n";

        oss << "\tend\n";
        oss << "end\n";
    }
    oss << "index_frac_bound = input('do you want to see frac bound and trace points? Input 1 for yes, and other any character for no.');\n";

    oss << "if (index_frac_bound == 1)\n";
    oss << "\t[m, ~] = size(Frac_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Frac_bond(:, 1), Frac_bond(:, 2), Frac_bond(:, 3), 'black', 'o', 'filled');\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "\t[m, ~] = size(Trace_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Trace_bond(:, 1), Trace_bond(:, 2), Trace_bond(:, 3), '*', 'LineWidth', 4);\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "end\n";

    oss << "\n\n";

    oss << "Top_bond = [];\n";
    oss << "Top_count = 1;\n";

    oss << "Bottom_bond = [];\n";
    oss << "Bottom_count = 1;\n";

    oss << "Front_bond = [];\n";
    oss << "Front_count = 1;\n";

    oss << "Back_bond = [];\n";
    oss << "Back_count = 1;\n";

    oss << "Left_bond = [];\n";
    oss << "Left_count = 1;\n";

    oss << "Right_bond = [];\n";
    oss << "Right_count = 1;\n";

    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        oss << "[m" << i << ", ~] = size(Frac_" << i << "_Pnt_attri_Mo_B);\n";
        oss << "for i = 1 : m" << i << "\n";
        oss << "\tfor j = 1 : 6\n";

        oss << "\t\tif (Frac_" << i << "_Pnt_attri_Mo_B(i,j) == 1 && j == 1)\n";
        oss << "\t\t\tTop_bond(Top_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tTop_count = Top_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n\n";

        oss << "\t\tif (Frac_" << i << "_Pnt_attri_Mo_B(i,j) == 1 && j == 2)\n";
        oss << "\t\t\tBottom_bond(Bottom_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tBottom_count = Bottom_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n\n";

        oss << "\t\tif (Frac_" << i << "_Pnt_attri_Mo_B(i,j) == 1 && j == 3)\n";
        oss << "\t\t\tFront_bond(Front_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tFront_count = Front_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n\n";

        oss << "\t\tif (Frac_" << i << "_Pnt_attri_Mo_B(i,j) == 1 && j == 4)\n";
        oss << "\t\t\tBack_bond(Back_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tBack_count = Back_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n\n";

        oss << "\t\tif (Frac_" << i << "_Pnt_attri_Mo_B(i,j) == 1 && j == 5)\n";
        oss << "\t\t\tLeft_bond(Left_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tLeft_count = Left_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n\n";

        oss << "\t\tif (Frac_" << i << "_Pnt_attri_Mo_B(i,j) == 1 && j == 6)\n";
        oss << "\t\t\tRight_bond(Right_count, :) = Frac_" << i << "_JXY3D(i, :);\n";
        oss << "\t\t\tRight_count = Right_count + 1;\n";
        oss << "\t\t\thold on;\n";
        oss << "\t\tend\n\n";

        oss << "\tend\n";
        oss << "end\n";
    }

    oss << "index_frac_bound = input('do you want to see model bound points? Input 1 for yes, and other any character for no.');\n";

    oss << "if (index_frac_bound == 1)\n";

    oss << "\t[m, ~] = size(Top_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Top_bond(:, 1), Top_bond(:, 2), Top_bond(:, 3), 'o', 'filled');\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "\t[m, ~] = size(Bottom_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Bottom_bond(:, 1), Bottom_bond(:, 2), Bottom_bond(:, 3), 'o', 'filled');\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "\t[m, ~] = size(Front_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Front_bond(:, 1), Front_bond(:, 2), Front_bond(:, 3), 'o', 'filled');\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "\t[m, ~] = size(Back_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Back_bond(:, 1), Back_bond(:, 2), Back_bond(:, 3), 'o', 'filled');\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "\t[m, ~] = size(Left_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Left_bond(:, 1), Left_bond(:, 2), Left_bond(:, 3), 'o', 'filled');\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "\t[m, ~] = size(Right_bond);\n";
    oss << "\tif (m > 0)\n";
    oss << "\t\tscatter3( Right_bond(:, 1), Right_bond(:, 2), Right_bond(:, 3), 'o', 'filled');\n";
    oss << "\t\thold on;\n";
    oss << "\tend\n";

    oss << "\tend\n";

    oss.close();
}

inline void Mesh_DFN::Ident_Model_Bound_Pnt(DFN::Domain dom)
{
    ///< dom.Model_domain;
    ///< Top-zmax, bottom-zmin, front-ymin, back-ymax, left-xmin, right-xmax
    this->Pnt_attri_Mo_B.resize(this->JXY_3D.size());
    for (size_t i = 0; i < this->JXY_3D.size(); ++i)
    {
        Pnt_attri_Mo_B[i].resize(this->JXY_3D[i].size());
        for (size_t j = 0; j < this->JXY_3D[i].size(); ++j)
        {
            double top_zmax = dom.Model_domain(0);
            double bottom_zmin = dom.Model_domain(1);
            double front_ymin = dom.Model_domain(2);
            double back_ymax = dom.Model_domain(3);
            double left_xmin = dom.Model_domain(4);
            double right_xmax = dom.Model_domain(5);

            Vector3d thisPnT = JXY_3D[i][j];
            if (abs(thisPnT(2) - top_zmax) < 0.001)
                Pnt_attri_Mo_B[i][j].If_model_top = true;

            if (abs(thisPnT(2) - bottom_zmin) < 0.001)
                Pnt_attri_Mo_B[i][j].If_model_bottom = true;

            if (abs(thisPnT(1) - front_ymin) < 0.001)
                Pnt_attri_Mo_B[i][j].If_model_front = true;

            if (abs(thisPnT(1) - back_ymax) < 0.001)
                Pnt_attri_Mo_B[i][j].If_model_back = true;

            if (abs(thisPnT(0) - left_xmin) < 0.001)
                Pnt_attri_Mo_B[i][j].If_model_left = true;

            if (abs(thisPnT(0) - right_xmax) < 0.001)
                Pnt_attri_Mo_B[i][j].If_model_right = true;
        }
    }
};

inline void Mesh_DFN::Find_repetitive_pnt(DFN::Domain dom)
{
    FEM_matrix_dimension = 0;
    Guide_frac.resize(JXY_3D.size());
    for (size_t i = 0; i < JXY_3D.size(); ++i)
    {
        Guide_frac[i].resize(JXY_3D[i].size());
        for (size_t j = 0; j < JXY_3D[i].size(); ++j)
        {
            Guide_frac[i][j] = Eigen::Vector3d::Zero();
        }
    }

    for (size_t i = 0; i < Trace_Tag.size(); ++i)
    {
        for (size_t j = 0; j < Trace_Tag[i].size(); ++j)
        {
            size_t PointID = Trace_Tag[i][j];

            if (Guide_frac[i][PointID](0) == 0)
            {
                FEM_matrix_dimension++;
                Vector3d thisPnt = JXY_3D[i][PointID];

                size_t FracID = Frac_Tag[i];

                for (std::set<size_t>::iterator it_k = dom.Fractures[FracID].Intersect_other_frac_after_trim.begin(); it_k != dom.Fractures[FracID].Intersect_other_frac_after_trim.end(); it_k++)
                {
                    for (size_t k = 0; k < Frac_Tag.size(); ++k)
                    {
                        if (*it_k == Frac_Tag[k])
                        {
                            for (size_t l = 0; l < Trace_Tag[k].size(); ++l)
                            {
                                size_t otherPntID = Trace_Tag[k][l];
                                Vector3d otherPnt = JXY_3D[k][otherPntID];
                                Vector3d KKO = thisPnt - otherPnt;
                                if (abs(KKO(0)) < 0.001 &&
                                    abs(KKO(1)) < 0.001 &&
                                    abs(KKO(2)) < 0.001)
                                {
                                    if (Guide_frac[k][otherPntID](0) == 0)
                                    {
                                        Guide_frac[k][otherPntID](0) = 1;       // repetitive point
                                        Guide_frac[k][otherPntID](1) = i;       // the i th frac
                                        Guide_frac[k][otherPntID](2) = PointID; // the PointID th point
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

}; // namespace DFN
