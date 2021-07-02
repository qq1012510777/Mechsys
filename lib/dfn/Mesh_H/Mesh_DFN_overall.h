#pragma once
#include "../DFN_H/Domain_WL.h"
#include "../Geometry_H/Point_3D.h"
#include "mat.h"

//
#include <gmsh.h>

typedef Eigen::Matrix<int, 1, 6> RowVector6i;

typedef struct Point_attribute
{
    bool If_model_top = false;
    bool If_model_bottom = false;
    bool If_model_front = false;
    bool If_model_back = false;
    bool If_model_left = false;
    bool If_model_right = false;
    //Top-zmax, bottom-zmin, front-ymin, back-ymax, left-xmin, right-xmax

    bool If_frac_bound = false;
    bool If_trace = false;
    //
    bool If_corner = false;
} P_attri;

namespace DFN
{
class Mesh_DFN_overall
{
public:
    std::vector<Vector3d> JXY_3D;
    //VectorXd PntTagLinear;
    vector<vector<Vector3d>> JXY_2D;

    std::vector<P_attri> Pnt_attri;

    std::vector<RowVector6i> JM;
    ///< each Vector3d is the ID values of
    ///< the three nodes

    std::vector<size_t> Frac_Tag;

    std::vector<std::vector<RowVector6i>> JM_Each_Frac;

    size_t NUM_of_NODES;
    size_t NUM_of_linear_NODES;

public:
    Mesh_DFN_overall();
    Mesh_DFN_overall(DFN::Domain dom, const double min_ele_edge, const double max_ele_edge);
    void Matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom);
    void Identify_point_attribute_and_2D_meshes(DFN::Domain dom);
    void Resort_pnts_orders();
    void Rotate_JXY_3D_to_2D(DFN::Domain dom);
};

inline Mesh_DFN_overall::Mesh_DFN_overall(){};

inline Mesh_DFN_overall::Mesh_DFN_overall(DFN::Domain dom, const double min_ele_edge, const double max_ele_edge)
{
    //---------------------------------------
#pragma omp critical
    {

        size_t VertsPntID = 1;
        size_t LineID = 1;
        size_t CurveLoopID = 1;
        size_t SurfaceID = 1;
        gmsh::initialize();
        gmsh::option::setNumber("General.Verbosity", 2); // default level is 5
        gmsh::model::add("t2");

        for (size_t i = 0; i < dom.Percolation_cluster[2].size(); ++i)
        {
            size_t ClusterID = dom.Percolation_cluster[2][i];

            for (size_t j = 0; j < dom.Listofclusters[ClusterID].size(); ++j)
            {
                size_t FracID = dom.Listofclusters[ClusterID][j];

                Frac_Tag.push_back(FracID);
                std::vector<int> Pointloop(dom.Fractures[FracID].Verts_trim.size());
                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                    Pointloop[k] = VertsPntID + k;

                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                {

                    gmsh::model::occ::addPoint(dom.Fractures[FracID].Verts_trim[k](0),
                                               dom.Fractures[FracID].Verts_trim[k](1),
                                               dom.Fractures[FracID].Verts_trim[k](2),
                                               0,
                                               VertsPntID);
                    VertsPntID++;
                }

                std::vector<int> curveloop(dom.Fractures[FracID].Verts_trim.size());
                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                    curveloop[k] = LineID + k;

                for (size_t k = 0; k < dom.Fractures[FracID].Verts_trim.size(); ++k)
                {
                    gmsh::model::occ::addLine(Pointloop[k],
                                              Pointloop[(k + 1) % dom.Fractures[FracID].Verts_trim.size()],
                                              LineID);
                    LineID++;
                }
                gmsh::model::occ::addCurveLoop(curveloop, CurveLoopID);

                std::vector<int> surfaceloop = {(int)CurveLoopID};
                gmsh::model::occ::addPlaneSurface(surfaceloop, SurfaceID);
                gmsh::model::occ::synchronize();

                CurveLoopID++;
                SurfaceID++;
            }
        }

        std::vector<std::pair<int, int>> input_entity(SurfaceID - 1);

        for (size_t i = 0; i < input_entity.size(); ++i)
            input_entity[i] = std::make_pair(2, i + 1);

        std::vector<std::pair<int, int>> out;
        std::vector<std::vector<std::pair<int, int>>> outmap;

        gmsh::model::occ::fragment(input_entity, input_entity, out, outmap);
        gmsh::model::occ::synchronize();

        gmsh::option::setNumber("Mesh.MeshSizeMin", min_ele_edge);
        gmsh::option::setNumber("Mesh.MeshSizeMax", max_ele_edge);

        gmsh::option::setNumber("Mesh.Algorithm", 5);

        //std::cout << "\033[31mstart meshing;\n\033[0m";
        gmsh::model::mesh::generate(2);
        //std::cout << "\033[31mfinish meshing;\n\033[0m";
        //--------
        double mw = 0;
        gmsh::option::getNumber("Mesh.NbNodes", mw);
        NUM_of_linear_NODES = mw;
        //-----------

        gmsh::option::setNumber("Mesh.ElementOrder", 2);
        gmsh::model::mesh::setOrder(2);

        std::vector<std::size_t> nodes;
        std::vector<double> coord, coordParam;
        gmsh::model::mesh::getNodes(nodes, coord, coordParam);
        NUM_of_NODES = coord.size() / 3;

        for (size_t i = 0; i < coord.size(); i += 3)
        {
            Vector3d A;
            A << coord[i], coord[i + 1], coord[i + 2];
            //cout << A.transpose() << endl;
            JXY_3D.push_back(A);
        }

        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 2, -1);

        for (size_t i = 0; i < elemNodeTags[0].size(); i += 6)
        {
            RowVector6i A;
            A << elemNodeTags[0][i] - 1,
                elemNodeTags[0][i + 3] - 1,
                elemNodeTags[0][i + 1] - 1,
                elemNodeTags[0][i + 4] - 1,
                elemNodeTags[0][i + 2] - 1,
                elemNodeTags[0][i + 5] - 1;
            //cout << A << endl;
            JM.push_back(A);
        }
        //gmsh::fltk::run();
        gmsh::clear();
        gmsh::finalize();

        this->Resort_pnts_orders();

        this->Identify_point_attribute_and_2D_meshes(dom);
        this->Rotate_JXY_3D_to_2D(dom);

    };
};

inline void Mesh_DFN_overall::Resort_pnts_orders()
{
    size_t PNT_NUM = this->JXY_3D.size();
    std::set<size_t> CornerID, MiddleID;
    for (size_t i = 0; i < this->JM.size(); ++i)
    {
        CornerID.insert(size_t(this->JM[i][0]));
        CornerID.insert(size_t(this->JM[i][2]));
        CornerID.insert(size_t(this->JM[i][4]));

        MiddleID.insert(size_t(this->JM[i][1]));
        MiddleID.insert(size_t(this->JM[i][3]));
        MiddleID.insert(size_t(this->JM[i][5]));
    }

    std::vector<size_t> old_vs_newID(PNT_NUM); // old is the index, new is the size_t element

    size_t _temp_ID_ = 0;
    for (std::set<size_t>::iterator its = CornerID.begin();
         its != CornerID.end(); ++its)
    {
        //its->second = _temp_ID_;
        old_vs_newID[*its] = _temp_ID_; // it* is the old ID
        _temp_ID_++;
    }

    for (std::set<size_t>::iterator its = MiddleID.begin();
         its != MiddleID.end(); ++its)
    {
        //its->second = _temp_ID_;
        old_vs_newID[*its] = _temp_ID_;
        _temp_ID_++;
    }

    std::vector<Vector3d> newJXY_3D(PNT_NUM);
    for (size_t i = 0; i < PNT_NUM; ++i)
    {
        size_t newID = old_vs_newID[i];
        newJXY_3D[newID] = this->JXY_3D[i];
    }

    std::vector<RowVector6i> newJM(this->JM.size());
    for (size_t i = 0; i < this->JM.size(); ++i)
    {
        for (size_t j = 0; j < 6; ++j)
        {
            size_t oldID = this->JM[i][j];
            size_t newID = old_vs_newID[oldID];
            newJM[i][j] = newID;
        }
    }

    this->JXY_3D = newJXY_3D;
    this->JM = newJM;
};

inline void Mesh_DFN_overall::Identify_point_attribute_and_2D_meshes(DFN::Domain dom)
{
    JM_Each_Frac.resize(Frac_Tag.size());
    //
    VectorXd Idx_u;
    Idx_u = Eigen::VectorXd::Zero(this->JXY_3D.size());

    //size_t PntTag_p = 0;

    //PntTagLinear = Eigen::VectorXd::Zero(this->JXY_3D.size());
    //PntTagLinear = PntTagLinear.array() - 1;

    this->Pnt_attri.resize(this->JXY_3D.size());

    for (size_t i = 0; i < this->JM.size(); ++i)
    {
        for (size_t j = 0; j < 6; ++j)
        {
            size_t PntID = this->JM[i](j);

            if (Idx_u[PntID] == 0)
            {
                Idx_u[PntID] = 1;
                Vector3d thisPnT = JXY_3D[PntID];
                DFN::Point_3D ThisPnt{thisPnT};

                // model bounds
                double top_zmax = dom.Model_domain(0);
                double bottom_zmin = dom.Model_domain(1);
                double front_ymin = dom.Model_domain(2);
                double back_ymax = dom.Model_domain(3);
                double left_xmin = dom.Model_domain(4);
                double right_xmax = dom.Model_domain(5);

                if (abs(thisPnT(2) - top_zmax) < 0.001)
                {
                    Pnt_attri[PntID].If_model_top = true;
                }

                if (abs(thisPnT(2) - bottom_zmin) < 0.001)
                    Pnt_attri[PntID].If_model_bottom = true;

                if (abs(thisPnT(1) - front_ymin) < 0.001)
                    Pnt_attri[PntID].If_model_front = true;

                if (abs(thisPnT(1) - back_ymax) < 0.001)
                    Pnt_attri[PntID].If_model_back = true;

                if (abs(thisPnT(0) - left_xmin) < 0.001)
                    Pnt_attri[PntID].If_model_left = true;

                if (abs(thisPnT(0) - right_xmax) < 0.001)
                    Pnt_attri[PntID].If_model_right = true;

                // frac bounds
                for (size_t k = 0; k < this->Frac_Tag.size(); ++k)
                {
                    size_t FracID_ = Frac_Tag[k];
                    if (ThisPnt.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true)
                    {
                        Pnt_attri[PntID].If_frac_bound = true;
                        break;
                    };
                }

                // frac trace
                for (std::map<std::pair<size_t, size_t>, std::pair<Vector3d, Vector3d>>::iterator its = dom.Intersections.begin();
                     its != dom.Intersections.end();
                     its++)
                {
                    std::vector<Vector3d> Trace_yu = {its->second.first, its->second.second};
                    if (ThisPnt.If_lies_on_a_line_seg(Trace_yu) == true)
                    {
                        Pnt_attri[PntID].If_trace = true;
                        break;
                    }
                }

                //corner
                if (j % 2 == 0)
                {
                    Pnt_attri[PntID].If_corner = true;
                    //PntTagLinear[PntID] = PntTag_p;
                    //PntTag_p++;
                }
            }
        }

        // which frac does this element lie on
        DFN::Point_3D ElePnt_1{JXY_3D[this->JM[i](0)]};
        DFN::Point_3D ElePnt_2{JXY_3D[this->JM[i](2)]};
        DFN::Point_3D ElePnt_3{JXY_3D[this->JM[i](4)]};

        for (size_t k = 0; k < this->Frac_Tag.size(); ++k)
        {
            size_t FracID_ = this->Frac_Tag[k];
            if ((ElePnt_1.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true || ElePnt_1.If_lies_within_a_polygon_3D(dom.Fractures[FracID_].Verts_trim) == true) &&
                (ElePnt_2.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true || ElePnt_2.If_lies_within_a_polygon_3D(dom.Fractures[FracID_].Verts_trim) == true) &&
                (ElePnt_3.If_lies_on_the_bounds_of_polygon(dom.Fractures[FracID_].Verts_trim) == true || ElePnt_3.If_lies_within_a_polygon_3D(dom.Fractures[FracID_].Verts_trim) == true))
            {
                JM_Each_Frac[k].push_back(this->JM[i]); //
                break;
            }

            if (k == this->Frac_Tag.size() - 1)
            {
                throw Error_throw_ignore("Cannot find which frac does this element lie on. In class Mesh_DFN_overall");
            }
        }
    }
};

void Mesh_DFN_overall::Matlab_plot(string FileKey_mat, string FileKey_m, DFN::Domain dom)
{
    const char *filename = FileKey_mat.c_str();
    MATFile *pMatFile;
    pMatFile = matOpen(filename, "w");

    if (!pMatFile)
    {
        throw Error_throw_ignore("cannot create mat file in class Mesh_DFN_overall\n");
    }

    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        //cout << "i: " << i << endl;
        size_t len = dom.Fractures[Frac_Tag[i]].Verts_trim.size(); // number of verts

        double *pData1;
        double *pData2;
        double *pData3;
        double *pData7;

        pData1 = (double *)mxCalloc(len, sizeof(double));
        pData2 = (double *)mxCalloc(len, sizeof(double));
        pData3 = (double *)mxCalloc(len, sizeof(double));
        pData7 = (double *)mxCalloc(this->JM_Each_Frac[i].size() * 6, sizeof(double));

        mxArray *pMxArray1;
        mxArray *pMxArray2;
        mxArray *pMxArray3;
        mxArray *pMxArray7;

        pMxArray1 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray2 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray3 = mxCreateDoubleMatrix(len, 1, mxREAL);
        pMxArray7 = mxCreateDoubleMatrix(this->JM_Each_Frac[i].size(), 6, mxREAL);

        if (!pMxArray1 || !pMxArray2 || !pMxArray3 || !pMxArray7)
        {
            throw Error_throw_ignore("cannot create pMxArray in class Mesh_DFN\n");
        }

        if (!pData1 || !pData2 || !pData3 || !pData7)
        {
            throw Error_throw_ignore("cannot create pData in class Mesh_DFN\n");
        }

        for (size_t j = 0; j < len; j++)
        {
            pData1[j] = dom.Fractures[Frac_Tag[i]].Verts_trim[j](0);
            pData2[j] = dom.Fractures[Frac_Tag[i]].Verts_trim[j](1);
            pData3[j] = dom.Fractures[Frac_Tag[i]].Verts_trim[j](2);
        }

        for (size_t j = 0; j < this->JM_Each_Frac[i].size() * 6; ++j)
        {
            size_t k, l;
            k = ceil(j / this->JM_Each_Frac[i].size()); // column
            l = j % this->JM_Each_Frac[i].size();       // row

            pData7[j] = this->JM_Each_Frac[i][l](k) + 1;
        }

        mxSetData(pMxArray1, pData1);
        mxSetData(pMxArray2, pData2);
        mxSetData(pMxArray3, pData3);
        mxSetData(pMxArray7, pData7);

        string ft = to_string(i);

        string Fracx = "Frac_" + ft + "_x";
        string Fracy = "Frac_" + ft + "_y";
        string Fracz = "Frac_" + ft + "_z";
        string JM_each = "JM_each_" + ft;

        const char *Fracx_s = Fracx.c_str();
        const char *Fracy_s = Fracy.c_str();
        const char *Fracz_s = Fracz.c_str();
        const char *JM_each_s = JM_each.c_str();

        matPutVariable(pMatFile, Fracx_s, pMxArray1);
        matPutVariable(pMatFile, Fracy_s, pMxArray2);
        matPutVariable(pMatFile, Fracz_s, pMxArray3);
        matPutVariable(pMatFile, JM_each_s, pMxArray7);

        mxFree(pData1);
        mxFree(pData2);
        mxFree(pData3);
        mxFree(pData7);
    }

    double *pData4;
    double *pData5;
    double *pData6;
    pData4 = (double *)mxCalloc(this->JXY_3D.size() * 3, sizeof(double));
    pData5 = (double *)mxCalloc(this->JM.size() * 6, sizeof(double));
    pData6 = (double *)mxCalloc(this->JXY_3D.size() * 9, sizeof(double));

    mxArray *pMxArray4;
    mxArray *pMxArray5;
    mxArray *pMxArray6;

    pMxArray4 = mxCreateDoubleMatrix(this->JXY_3D.size(), 3, mxREAL);
    pMxArray5 = mxCreateDoubleMatrix(this->JM.size(), 6, mxREAL);
    pMxArray6 = mxCreateDoubleMatrix(this->JXY_3D.size(), 9, mxREAL);

    if (!pMxArray4 || !pMxArray5 || !pMxArray6)
    {
        throw Error_throw_ignore("cannot create pMxArray in class Mesh_DFN\n");
    }

    if (!pData4 || !pData5 || !pData6)
    {
        throw Error_throw_ignore("cannot create pData in class Mesh_DFN\n");
    }

    for (size_t j = 0; j < this->JXY_3D.size() * 3; ++j)
    {
        size_t k, l;
        k = ceil(j / this->JXY_3D.size()); // column
        l = j % this->JXY_3D.size();       // row

        pData4[j] = this->JXY_3D[l](k);
        //cout << pData4[j] << endl;
    }

    for (size_t j = 0; j < this->JM.size() * 6; ++j)
    {
        size_t k, l;
        k = ceil(j / this->JM.size()); // column
        l = j % this->JM.size();       // row

        pData5[j] = this->JM[l](k) + 1;
    }

    for (size_t j = 0; j < this->JXY_3D.size() * 9; ++j)
    {
        size_t k, l;
        k = ceil(j / this->JXY_3D.size()); // column
        l = j % this->JXY_3D.size();       // row

        if (k == 0)
            pData6[j] = (int)this->Pnt_attri[l].If_model_top;
        else if (k == 1)
            pData6[j] = (int)this->Pnt_attri[l].If_model_bottom;
        else if (k == 2)
            pData6[j] = (int)this->Pnt_attri[l].If_model_front;
        else if (k == 3)
            pData6[j] = (int)this->Pnt_attri[l].If_model_back;
        else if (k == 4)
            pData6[j] = (int)this->Pnt_attri[l].If_model_left;
        else if (k == 5)
            pData6[j] = (int)this->Pnt_attri[l].If_model_right;
        else if (k == 6)
            pData6[j] = (int)this->Pnt_attri[l].If_frac_bound;
        else if (k == 7)
            pData6[j] = (int)this->Pnt_attri[l].If_trace;
        else if (k == 8)
            pData6[j] = (int)this->Pnt_attri[l].If_corner;
        //cout << pData4[j] << endl;
    }

    mxSetData(pMxArray4, pData4);
    mxSetData(pMxArray5, pData5);
    mxSetData(pMxArray6, pData6);

    const char *FracJXY3D_s = "Frac_JXY3D";
    const char *FracJM_s = "Frac_JM";
    const char *PntAttri_s = "PntAttri_s";

    matPutVariable(pMatFile, FracJXY3D_s, pMxArray4);
    matPutVariable(pMatFile, FracJM_s, pMxArray5);
    matPutVariable(pMatFile, PntAttri_s, pMxArray6);

    mxFree(pData4);
    mxFree(pData5);
    mxFree(pData6);

    matClose(pMatFile);

    // m file
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;";
    oss << "%%---matrix dimension is " << NUM_of_NODES * 2 + NUM_of_linear_NODES << "--%%\n"
        << endl;
    oss << "load('" << FileKey_mat << "');\n";

    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        oss << "fill3([Frac_" << i << "_x; Frac_" << i << "_x(1,1)], [Frac_" << i << "_y; Frac_" << i << "_y(1,1)], [Frac_" << i << "_z; Frac_" << i << "_z(1,1)], [rand rand rand]);\nhold on;\n";
    }

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

    oss << "[m, ~] = size(Frac_JXY3D);\n";
    oss << "P = patch('Vertices', Frac_JXY3D, 'Faces', Frac_JM, 'FaceVertexCData', zeros(m, 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);\n";
    oss << "hold on;\n";

    oss << "index_frac_bound = input('do you want to see special points? Input 1 for yes, and other any character for no.');\n";

    oss << "if (index_frac_bound == 1)\n";
    oss << "\t[m, ~] = size(" << FracJXY3D_s << ");\n";

    oss << "\t\tTOP = [];\n";
    oss << "\t\tBOTTOM = [];\n";
    oss << "\t\tFRONT = [];\n";
    oss << "\t\tBACK = [];\n";
    oss << "\t\tLEFT = [];\n";
    oss << "\t\tRIGHT = [];\n";
    oss << "\t\tBOUND = [];\n";
    oss << "\t\tTRACE = [];\n";
    oss << "\t\tCORNER = [];\n";

    oss << "\tfor i = 1:m;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 1) == 1)\n";
    oss << "\t\t\tTOP = [TOP i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 2) == 1)\n";
    oss << "\t\t\tBOTTOM = [BOTTOM i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 3) == 1 & " << PntAttri_s << "(i, 1) ~= 1 & " << PntAttri_s << "(i, 2) ~= 1)\n";
    oss << "\t\t\tFRONT= [FRONT i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 4) == 1 & " << PntAttri_s << "(i, 1) ~= 1 & " << PntAttri_s << "(i, 2) ~= 1)\n";
    oss << "\t\t\tBACK= [BACK i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 5) == 1 & " << PntAttri_s << "(i, 1) ~= 1 & " << PntAttri_s << "(i, 2) ~= 1)\n";
    oss << "\t\t\tLEFT= [LEFT i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 6) == 1 & " << PntAttri_s << "(i, 1) ~= 1 & " << PntAttri_s << "(i, 2) ~= 1)\n";
    oss << "\t\t\tRIGHT= [RIGHT i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 7) == 1 & " << PntAttri_s << "(i, 1) ~= 1 & " << PntAttri_s << "(i, 2) ~= 1)\n";
    oss << "\t\t\tBOUND= [BOUND i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 8) == 1)\n";
    oss << "\t\t\tTRACE= [TRACE i];\n";
    oss << "\t\tend;\n";

    oss << "\t\tif (" << PntAttri_s << "(i, 9) == 1)\n";
    oss << "\t\t\tCORNER= [CORNER i];\n";
    oss << "\t\tend;\n";

    oss << "\tend;\n";
    oss << "end;\n";

    oss << "if (index_frac_bound == 1)\n";
    ;
    oss << "\tindex_frac_bound2 = input('do you want to see top and bottom points? Input 1 for yes, and other any character for no.');\n";

    oss << "\tif (index_frac_bound2 == 1)\n";
    oss << "\t\tscatter3(Frac_JXY3D([TOP BOTTOM], 1), Frac_JXY3D([TOP BOTTOM], 2), Frac_JXY3D([TOP BOTTOM], 3), 'o', 'green', 'filled');\n";
    oss << "\tend;\n";

    oss << "\tindex_frac_bound2 = input('do you want to see lateral surface points? Input 1 for yes, and other any character for no.');\n";

    oss << "\tif (index_frac_bound2 == 1)\n";
    oss << "\t\tscatter3(Frac_JXY3D([FRONT BACK LEFT RIGHT], 1), Frac_JXY3D([FRONT BACK LEFT RIGHT], 2), Frac_JXY3D([FRONT BACK LEFT RIGHT], 3), 'o', 'blue', 'filled');\n";
    oss << "\tend;\n";

    oss << "\tindex_frac_bound2 = input('do you want to see frac bound points? Input 1 for yes, and other any character for no.');\n";

    oss << "\tif (index_frac_bound2 == 1)\n";
    oss << "\t\tscatter3(Frac_JXY3D([BOUND], 1), Frac_JXY3D([BOUND], 2), Frac_JXY3D([BOUND], 3), 'o', 'red', 'filled');\n";
    oss << "\tend;\n";

    oss << "\tindex_frac_bound2 = input('do you want to see trace points? Input 1 for yes, and other any character for no.');\n";

    oss << "\tif (index_frac_bound2 == 1)\n";
    oss << "\t\tscatter3(Frac_JXY3D([TRACE], 1), Frac_JXY3D([TRACE], 2), Frac_JXY3D([TRACE], 3), 'o', 'cyan', 'filled');\n";
    oss << "\tend;\n";

    oss << "\tindex_frac_bound2 = input('do you want to see corner points? Input 1 for yes, and other any character for no.');\n";

    oss << "\tif (index_frac_bound2 == 1)\n";
    oss << "\t\tscatter3(Frac_JXY3D([CORNER], 1), Frac_JXY3D([CORNER], 2), Frac_JXY3D([CORNER], 3), 'o', 'black', 'filled');\n";
    oss << "\tend;\n";
    oss << "end;\n";

    oss << "index_frac_bound3 = input('do you want to check topo structure for each frac? Input 1 for yes, and other any character for no.');\n";

    oss << "if (index_frac_bound3 == 1)\n";
    oss << "\tfigure(2);\n";
    oss << "\ttitle('topo check');\n";
    oss << "\tview(3);\n";
    oss << "\t[m, ~] = size(Frac_JXY3D);\n";
    for (size_t i = 0; i < this->Frac_Tag.size(); ++i)
    {
        oss << "\tP = patch('Vertices', Frac_JXY3D, 'Faces', JM_each_" << i << ", 'FaceVertexCData', ones(m, 1) .* " << 2 * i << ", 'FaceColor', 'interp', 'EdgeAlpha', 0.9);\n";
        oss << "\thold on;\n";
    }
    oss << "end;\n";

    oss << "\n\n\nfigure(1)\n"
        << endl;
    oss << "[m, ~] = size(Frac_JXY3D);\n";
    oss << "for i = 1:m\n";
    oss << "\ttext(Frac_JXY3D(i, 1), Frac_JXY3D(i, 2), Frac_JXY3D(i, 3), num2str(i), 'FontSize', 11);\n";
    oss << "end\n";

    //oss << "\n\n\nfigure(3)\n"
    //  << endl;
    /*
    for (size_t i = 0; i < Frac_Tag.size(); ++i)
    {
        oss << "fill3([Frac_" << i << "_x; Frac_" << i << "_x(1,1)], [Frac_" << i << "_y; Frac_" << i << "_y(1,1)], [Frac_" << i << "_z; Frac_" << i << "_z(1,1)], [rand rand rand]);\nhold on;\n";
    }
    for (int i = 0; i < PntTagLinear.size(); ++i)
    {
        if (PntTagLinear[i] != -1)
        {
            oss << "text(" << JXY_3D[i](0) << ", " << JXY_3D[i](1) << ", " << JXY_3D[i](2) << ", num2str(" << PntTagLinear[i] + 1 << "), 'FontSize', 11);\nhold on;\n";
        }
    }
    */

    oss.close();
};

inline void Mesh_DFN_overall::Rotate_JXY_3D_to_2D(DFN::Domain dom)
{
    JXY_2D.resize(JM_Each_Frac.size());

    for (size_t i = 0; i < JM_Each_Frac.size(); ++i)
    {
        JXY_2D[i].resize(this->JXY_3D.size());

        VectorXd PNT_INDICATOR = Eigen::VectorXd::Zero(this->JXY_3D.size());

        size_t Frac_Tag = this->Frac_Tag[i];

        DFN::Polygon_convex_3D poly{dom.Fractures[Frac_Tag].Verts_trim};
        std::vector<Vector3d> verts1;
        DFN::Rotate_to_horizontal R1{poly.Corners, verts1};

        for (size_t j = 0; j < JM_Each_Frac[i].size(); ++j)
        {
            for (size_t k = 0; k < 6; ++k)
            {
                size_t PNT_ID = JM_Each_Frac[i][j][k];
                if (PNT_INDICATOR[PNT_ID] == 0)
                {
                    PNT_INDICATOR[PNT_ID] = 1;

                    bool UY;
                    std::vector<Vector3d> jxy_3d(1), verts2;
                    jxy_3d[0] = this->JXY_3D[PNT_ID];
                    R1.Rotate_other_pnts(jxy_3d, verts2, UY);
                    if (UY == false)
                    {
                        cout << "Rotate pnts to 2D failed! In class 'MESH_DFN_overall'\n";
                        throw Error_throw_ignore("Rotate pnts to 2D failed! In class 'MESH_DFN_overall'!\n");
                    }

                    JXY_2D[i][PNT_ID] = verts2[0];
                }
            }
        }
    }
};

}; // namespace DFN
