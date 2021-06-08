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
    std::vector<size_t> Num_pressure;

    std::vector<std::vector<RowVector6i>> JM;
    ///< each Vector3d is the ID values of
    ///< the three nodes

    std::vector<vector<bool>> IfCorner;

    std::vector<std::vector<size_t>> Trace_Tag; // record point Tags belonging to traces of each fracture
    std::vector<std::vector<Eigen::Vector4d>> Guide_frac;
    ///< (0): an index referring if this is
    //a repetitive point;
    //(1): i; (2): j;(if repetitive)

    size_t NUM_of_NODES;
    size_t NUM_of_linear_NODES;

    std::map<std::pair<size_t, size_t>, std::vector<Vector3d>> Trace_been_split;

public:
    Mesh_DFN(DFN::Domain dom, double subtrace, double subpolygon);
    void Create_structure_of_a_frac(DFN::Domain dom, const size_t FracID, std::vector<Vector3d> &Verts, vector<std::pair<Vector3d, Vector3d>> &Traces, DFN::Remove_unnecess_frac Remove_f, std::vector<std::pair<size_t, size_t>> &Trace_ID);
    void Matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom);
    void Ident_Model_Bound_Pnt(DFN::Domain dom);
    void Find_repetitive_pnt(DFN::Domain dom);
    void Check_overall_pnt_sets(DFN::Domain dom);
    bool If_this_fracture_intersects_top_and_or_bot(std::vector<Vector3d> Verts_frac, double &toplength, double &botlength, DFN::Domain dom);
    void Identify_splitting_traces(const std::vector<std::pair<size_t, size_t>> Trace_ID, const std::vector<Vector3d> Traces_frac_f_rotated,
                                   std::map<std::pair<size_t, size_t>, std::vector<Vector2d>> &Trace_been_split_2D);
};

inline Mesh_DFN::Mesh_DFN(DFN::Domain dom, double subtrace, double subpolygon)
{
    double toplength = 0;
    double botlength = 0;
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
                std::vector<std::pair<size_t, size_t>> Trace_ID;
                this->Create_structure_of_a_frac(dom, FracID, Verts_frac, Traces_frac, Remove_f, Trace_ID);

                this->If_this_fracture_intersects_top_and_or_bot(Verts_frac, toplength, botlength, dom);

                DFN::Polygon_convex_3D poly{Verts_frac};
                poly.Optimize();
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
                //std::vector<DFN::Line_seg_2D> Traces_Line;

                for (size_t k = 0; k < Traces_frac_f_rotated.size(); k = k + 2)
                {
                    size_t pd = k / 2;
                    Traces_Line[pd].Re_constructor(Vector2d{Traces_frac_f_rotated[k](0), Traces_frac_f_rotated[k](1)}, Vector2d{Traces_frac_f_rotated[k + 1](0), Traces_frac_f_rotated[k + 1](1)});
                }

                ///--------------remove the trace that has been split and add the split trace segments
                /*
                for (size_t k = 0; k < Trace_ID.size(); ++k)
                {
                    std::map<std::pair<size_t, size_t>, std::vector<Vector3d>>::iterator ity = Trace_been_split.find(Trace_ID[k]);

                    if (ity == Trace_been_split.end())
                    {
                        DFN::Line_seg_2D Trace_ER{
                            Vector2d{Traces_frac_f_rotated[k * 2](0), Traces_frac_f_rotated[k * 2](1)},
                            Vector2d{Traces_frac_f_rotated[k * 2 + 1](0), Traces_frac_f_rotated[k * 2 + 1](1)},
                        };
                        Traces_Line.push_back(Trace_ER);
                    }
                    else
                    {
                        //
                        std::vector<Vector3d> Trace_SEG_3D = ity->second;
                        std::vector<Vector3d> Trace_SEG_2D = ity->second;
                        //rotate
                        if (poly.Beta > 0.0001)
                        {
                            DFN::Rotation_verts R3(Trace_SEG_3D, R_angle_temp1, Q_axis_1, temp1, Trace_SEG_2D);
                        }

                        //
                        for (size_t l = 0; l < Trace_SEG_2D.size() - 1; ++l)
                        {
                            DFN::Line_seg_2D Trace_ER{
                                Vector2d{Trace_SEG_2D[l](0), Trace_SEG_2D[l](1)},
                                Vector2d{Trace_SEG_2D[l + 1](0), Trace_SEG_2D[l + 1](1)},
                            };
                            Traces_Line.push_back(Trace_ER);
                        }
                    }
                }
                */
                //-------------------------------------------------------------------------------------

                // there should be a pre_handle code to divide polygon2D if extream traces exist
                // extream means: the trace nearly intersect edges but actually not
                // and very tiny trace but not a point

                DFN::Polygon_convex_2D polyframe{Verts_frac_rotated};

                DFN::Polygon_convex_2D_with_traces poly2D{polyframe, Traces_Line};

                DFN::Splitting_Polygon_convex_2D_with_traces splitting_poly2D{poly2D, subtrace, Traces_Line};
                //splitting_poly2D.Matlab_plot("D333frac.m");
                std::vector<std::vector<std::pair<size_t, Eigen::Vector2d>>> neigh_shared_A;

                size_t NO_Nodes_p_A = 0;

                DFN::Mesh_Polygon_2D_with_traces mesh{
                    splitting_poly2D.Pnt_sets,
                    splitting_poly2D.Topo_struct,
                    splitting_poly2D.Topo_struct_attri,
                    subpolygon,
                    neigh_shared_A,
                    NO_Nodes_p_A,
                    polyframe,
                    splitting_poly2D,
                    Traces_Line};

                Num_pressure.push_back(NO_Nodes_p_A);
                this->Pnt_attribute.push_back(mesh.Pnt_attribute);

                this->JXY.push_back(mesh.JXY);
                this->JM.push_back(mesh.JM);
                this->Trace_Tag.push_back(mesh.Trace_Tag);
                this->IfCorner.push_back(mesh.IfCorner);

                //-----------------
                // need a code to identify split traces
                std::map<std::pair<size_t, size_t>, std::vector<Vector2d>> Trace_been_split_2D;
                this->Identify_splitting_traces(Trace_ID, Traces_frac_f_rotated, Trace_been_split_2D);
                std::vector<std::vector<Vector3d>> TMP_trace_splitting_2D(Trace_been_split_2D.size()), TMP_trace_splitting_3D(Trace_been_split_2D.size());

                size_t k_idx = 0;
                for (std::map<std::pair<size_t, size_t>, std::vector<Vector2d>>::iterator iters = Trace_been_split_2D.begin();
                     iters != Trace_been_split_2D.end(); ++iters)
                {
                    TMP_trace_splitting_2D[k_idx].resize(iters->second.size());
                    TMP_trace_splitting_3D[k_idx].resize(iters->second.size());

                    for (size_t k = 0; k < iters->second.size(); ++k)
                    {
                        TMP_trace_splitting_2D[k_idx][k] << iters->second[k](0), iters->second[k](1), 0;
                        TMP_trace_splitting_3D[k_idx][k] << iters->second[k](0), iters->second[k](1), 0;
                    }

                    k_idx++;
                }
                //----------------

                std::vector<Vector3d> Verts_2D_tmp(mesh.JXY.size());
                for (size_t k = 0; k < Verts_2D_tmp.size(); ++k)
                    Verts_2D_tmp[k] << mesh.JXY[k](0), mesh.JXY[k](1), z_coord; // move back
                std::vector<Vector3d> Verts_rotated_back(mesh.JXY.size());

                if (poly.Beta > 0.0001)
                {
                    DFN::Rotation_verts R2(Verts_2D_tmp, -R_angle_temp1, Q_axis_1, temp1, Verts_rotated_back);

                    for (size_t k = 0; k < TMP_trace_splitting_2D.size(); ++k)
                        DFN::Rotation_verts R3(TMP_trace_splitting_2D[k], -R_angle_temp1, Q_axis_1, temp1, TMP_trace_splitting_3D[k]);
                }
                else
                {
                    Verts_rotated_back = Verts_2D_tmp;
                }
                this->JXY_3D.push_back(Verts_rotated_back);

                k_idx = 0;
                for (std::map<std::pair<size_t, size_t>, std::vector<Vector2d>>::iterator iters = Trace_been_split_2D.begin();
                     iters != Trace_been_split_2D.end(); ++iters)
                {
                    std::pair<std::pair<size_t, size_t>, std::vector<Vector3d>> tmp_pair;
                    tmp_pair.first = iters->first;
                    tmp_pair.second = TMP_trace_splitting_3D[k_idx];
                    this->Trace_been_split.insert(tmp_pair);

                    k_idx++;
                }
                //-----------------------------
            }
        }
    }

    this->Ident_Model_Bound_Pnt(dom);
    this->Find_repetitive_pnt(dom);
    //cout << "toplength: " << toplength << endl;
    //cout << "botlength: " << botlength << endl;

    if (toplength / subpolygon < 5)
    {
        throw Error_throw_ignore("Inlet length is too small!\n");
    }

    if (botlength / subpolygon < 5)
    {
        throw Error_throw_ignore("Outlet length is too small!\n");
    }
};

inline void Mesh_DFN::Create_structure_of_a_frac(DFN::Domain dom,
                                                 const size_t FracID,
                                                 std::vector<Vector3d> &Verts,
                                                 std::vector<std::pair<Vector3d, Vector3d>> &Traces,
                                                 DFN::Remove_unnecess_frac Remove_f,
                                                 std::vector<std::pair<size_t, size_t>> &Trace_ID)
{
    //

    // fracture vertices
    Verts.resize(dom.Fractures[FracID].Verts_trim.size());
    for (size_t i = 0; i < Verts.size(); ++i)
        Verts[i] = dom.Fractures[FracID].Verts_trim[i];

    // traces
    for (size_t i = 0; i < Remove_f.Listofclusters_remove_fracs.size(); ++i)
    {
        for (size_t j = 0; j < Remove_f.Listofclusters_remove_fracs[i].size(); ++j)
        {
            if ((size_t)Remove_f.Listofclusters_remove_fracs[i][j] != FracID)
            {
                size_t other_fracID = Remove_f.Listofclusters_remove_fracs[i][j];
                size_t ID_1 = FracID < other_fracID ? FracID : other_fracID;
                size_t ID_2 = FracID > other_fracID ? FracID : other_fracID;

                std::pair<size_t, size_t> Trace_key = std::make_pair(ID_1, ID_2);
                std::map<std::pair<size_t, size_t>, std::pair<Vector3d, Vector3d>>::iterator iter = dom.Intersections.find(Trace_key);

                if (iter != dom.Intersections.end())
                {
                    Traces.push_back(std::make_pair(dom.Intersections[Trace_key].first, dom.Intersections[Trace_key].second));
                    Trace_ID.push_back(std::make_pair(ID_1, ID_2));
                }
            }
        }
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
        throw Error_throw_ignore("cannot create mat file in class Mesh_DFN\n");
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
            throw Error_throw_ignore("cannot create pMxArray in class Mesh_DFN\n");
        }

        if (!pData1 || !pData2 || !pData3 || !pData4 || !pData5 || !pData6 || !pData7)
        {
            throw Error_throw_ignore("cannot create pData in class Mesh_DFN\n");
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
    oss << "%% matrix dimension is " << this->NUM_of_NODES * 2 + this->NUM_of_linear_NODES << endl;
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
    NUM_of_NODES = 0;
    NUM_of_linear_NODES = 0;

    Guide_frac.resize(JXY_3D.size());
    for (size_t i = 0; i < JXY_3D.size(); ++i)
    {
        Guide_frac[i].resize(JXY_3D[i].size());
        for (size_t j = 0; j < JXY_3D[i].size(); ++j)
        {
            Guide_frac[i][j] = Eigen::Vector4d::Zero();
        }
    }

    for (size_t i = 0; i < JXY.size(); ++i)
    {
        for (size_t j = 0; j < JXY[i].size(); ++j)
        {
            size_t PointID = j;

            if (Guide_frac[i][PointID](0) == 0)
            {
                Guide_frac[i][PointID](0) = NUM_of_NODES;
                NUM_of_NODES++;

                if (this->IfCorner[i][PointID] == true)
                {
                    Guide_frac[i][PointID](3) = NUM_of_linear_NODES;
                    NUM_of_linear_NODES++;
                }
                else
                {
                    Guide_frac[i][PointID](3) = -2; // meaning that it is not a pressure point / corner
                }

                if (Pnt_attribute[i][j].If_trace == true) // if it is a trace point
                {

                    Vector3d thisPnt = JXY_3D[i][PointID];

                    size_t FracID = Frac_Tag[i];
                    size_t NumOfShareTime = 0;
                    for (std::set<size_t>::iterator it_k = dom.Fractures[FracID].Intersect_other_frac_after_trim.begin();
                         it_k != dom.Fractures[FracID].Intersect_other_frac_after_trim.end();
                         it_k++)
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
                                    if ((abs(KKO(0)) < 0.1 &&
                                         abs(KKO(1)) < 0.1 &&
                                         abs(KKO(2)) < 0.1) ||
                                        KKO.norm() < 0.1)
                                    {
                                        NumOfShareTime++;
                                        if (Guide_frac[k][otherPntID](0) == 0)
                                        {
                                            Guide_frac[k][otherPntID](0) = -1;      // repetitive point
                                            Guide_frac[k][otherPntID](1) = i;       // the i th frac
                                            Guide_frac[k][otherPntID](2) = PointID; // the PointID th point

                                            if (IfCorner[i][PointID] == true)
                                                Guide_frac[k][otherPntID](3) = -1;
                                            else
                                                Guide_frac[k][otherPntID](3) = -2; // not a corner point
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (NumOfShareTime == 0)
                    {
                        cout << "cannot find the shared point in another intersecting fracture, in class Mesh_DFN!\n";
                        //cout << "FracID: " << i << endl;
                        //cout << "This point (3D) is:\n";
                        //cout << "\t[" << thisPnt(0) << ", 18], [" << thisPnt(1) << ", 18], [" << thisPnt(2) << ", 18]" << endl;
                        //cout << "Output the mesh file (mesh_check_pnt.m) now to check!\n";
                        //this->Matlab_plot("mesh_check_pnt.mat", "mesh_check_pnt.m", dom);
                        //dom.PlotMatlab_DFN_and_Intersection("tdfn01_DFN_and_Intersections2_V.m");
                        //dom.PLotMatlab_DFN_Cluster_along_a_direction("tdfn01_DFN_Z_clusters.m", "z");
                        throw Error_throw_ignore("cannot find the shared point in another intersecting fracture, in class Mesh_DFN!\n");
                    }
                }
            }
        }
    }
}

inline bool Mesh_DFN::If_this_fracture_intersects_top_and_or_bot(std::vector<Vector3d> Verts_frac, double &toplength, double &botlength, DFN::Domain dom)
{
    bool intersect = false;
    for (size_t i = 0; i < Verts_frac.size(); ++i)
    {
        double z1 = Verts_frac[i](2), z2 = Verts_frac[(i + 1) % Verts_frac.size()](2);

        if (abs(z1 - dom.Model_domain(0)) < 0.001 && abs(z2 - dom.Model_domain(0)) < 0.001)
        {
            intersect = true;
            Vector3d asd = Verts_frac[i] - Verts_frac[(i + 1) % Verts_frac.size()];

            toplength += asd.norm();
        }

        if (abs(z1 - dom.Model_domain(1)) < 0.001 && abs(z2 - dom.Model_domain(1)) < 0.001)
        {
            intersect = true;
            Vector3d asd = Verts_frac[i] - Verts_frac[(i + 1) % Verts_frac.size()];

            botlength += asd.norm();
        }
    }

    return intersect;
};

inline void Mesh_DFN::Identify_splitting_traces(const std::vector<std::pair<size_t, size_t>> Trace_ID, const std::vector<Vector3d> Traces_frac_f_rotated,
                                                std::map<std::pair<size_t, size_t>, std::vector<Vector2d>> &Trace_been_split_2D)
{
    //-----------
    for (size_t i = 0; i < Trace_ID.size(); ++i)
    {
        std::map<std::pair<size_t, size_t>, std::vector<Vector3d>>::iterator iter = this->Trace_been_split.find(Trace_ID[i]);

        if (iter == this->Trace_been_split.end()) // means that the splitting trace segments are not been recorded yet
        {
            std::pair<std::pair<size_t, size_t>, std::vector<Vector2d>> tmp_trace_seg;
            tmp_trace_seg.first = Trace_ID[i];

            DFN::Line_seg_2D Thisline{Vector2d{Traces_frac_f_rotated[i * 2](0), Traces_frac_f_rotated[i * 2](1)},
                                      Vector2d{Traces_frac_f_rotated[i * 2 + 1](0), Traces_frac_f_rotated[i * 2 + 1](1)}};

            for (size_t j = 0; j < this->JXY[this->JXY.size() - 1].size(); ++j)
            {
                DFN::Point_2D Thepnt{this->JXY[this->JXY.size() - 1][j]};
                bool AKI = Thepnt.If_lies_on_a_line_seg(Thisline.Point);
                if (AKI == true)
                {
                    tmp_trace_seg.second.push_back(Thepnt.Coordinate);
                }
            }

            //--------minus the first end
            for (size_t j = 0; j < tmp_trace_seg.second.size(); ++j)
                tmp_trace_seg.second[j] -= Vector2d{Traces_frac_f_rotated[i * 2](0), Traces_frac_f_rotated[i * 2](1)};

            //sort
            sort(tmp_trace_seg.second.begin(), tmp_trace_seg.second.end(), compVector2d);

            // plus back
            for (size_t j = 0; j < tmp_trace_seg.second.size(); ++j)
                tmp_trace_seg.second[j] += Vector2d{Traces_frac_f_rotated[i * 2](0), Traces_frac_f_rotated[i * 2](1)};

            Trace_been_split_2D.insert(tmp_trace_seg);
        }
    }
};

}; // namespace DFN
