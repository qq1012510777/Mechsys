#pragma once
//#include "../Geometry_H/Intersection_Frac.h"
#include "../Geometry_H/Intersection_Frac_boost.h"
#include "../Geometry_H/Intersection_between_polygon_and_3D_box.h"
#include "../Graph_WL_H/Graph_WL.h"
#include "../HDF5_DFN/HDF5_DFN.h"
#include "../MATLAB_DATA_API/MATLAB_DATA_API.h"
#include "Fracture_WL.h"
#include "mat.h"
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace Eigen;

namespace DFN
{
class Domain
{
public:
    //Data
    std::vector<Fracture> Fractures;                                                  ///< std::vector of Fractures, this array stores all generated fractures
    std::vector<size_t> Connections;                                                  ///< std::vector of Fracture connectivity, record the tag / ID of fractures that are connected
    std::vector<std::vector<size_t>> Listofclusters;                                  ///< List of fractures per cluster
    std::vector<std::vector<size_t>> Percolation_cluster;                             ///< three orders, N dimensions; first order means x direction, second means y, third means z; alone a percolation direction, there might be zero or more clusters; each element is the subscript of array (Listofclusters)
    std::map<std::pair<size_t, size_t>, std::pair<Vector3d, Vector3d>> Intersections; ///< Map of line intersection between pairs of fractures
    Vector6d Model_domain;                                                            ///< Top-zmax, bottom-zmin, front-ymin, back-ymax, left-xmin, right-xmax

    std::vector<bool> Percolation_status = {false, false, false}; ///< if the model is percolative along x, y ,and z
    double n_I;                                                   ///< Average number of intersections per fracture
    double P30;
    double P32_total;

    Vector3d P30_connected; // linear density of percolating clusters
    Vector3d P32_connected; // areal density of percolating clusters

    Vector3d Ratio_of_P32; ///< areal probability of a fracture belonging to percolating clusters
    Vector3d Ratio_of_P30; ///< linear probability of a fracture belonging to percolating clusters

    double P30_largest_cluster; // linear density of largest clusters
    double P32_largest_cluster; // areal density of largest clusters

    double Last_frac_size;

    std::vector<Fracture> Surfaces;                                                     ///< model surface
    std::vector<size_t> Connections_S;                                                  ///< thos fractures intersect with surfaces
    std::map<std::pair<size_t, size_t>, std::pair<Vector3d, Vector3d>> Intersections_S; ///< Map of line intersection between fractures and surfaces

    bool mode_2D = false;

public:
    void Re_identify_intersection_considering_trimmed_frac();
    // method
    void Create_whole_model(const size_t n,
                            const std::vector<double> DenWeight,
                            gsl_rng *random_seed,
                            const Vector6d model_size,
                            const string str_ori,
                            const string str_frac_size,
                            const std::vector<Vector2d> array11,
                            std::vector<Vector4d> array12,
                            const std::vector<Vector7d> array13,
                            const string conductivity_distri);

    void Model_set(const Vector6d model_size);
    ///< define model domain

    void AddSquareFracture(size_t Tag, Fracture &c);
    ///< Adding a square fracture

    bool Intersect_A(const Fracture F1, const Fracture F2);
    ///< JUST Identify if fractures
    // are connected, but no intersection points
    // are returned

    void Modify_trimmed_fractures_attribute(Fracture &F2);

    void If_fracture_intersect_boundary(Fracture &F2);
    // if a fracture intersects boundary, then
    // label it

    bool Intersect(Fracture &F1, Fracture &F2, const bool if_use_trim_frac);
    ///< Function to check if two fractures
    // intersect, and return intersection

    void Clusters();
    ///< Check the connectivity array
    // to form the clusters

    void Determine_max_cluster();
    //

    void Average_number_of_intersections_per_fracture();
    ///< Average_number_of_intersections_per_fracture

    void Determine_excluded_volume(const string str_ori, const string str_frac_size, double alpha_g = 0, double kappa = 0, double mean_i = 0, double var_i = 0, double min_R_i = 0, double max_R_i = 0);
    // do not use

    void Identify_percolation_clusters();
    //

    void Connectivity_analysis(); ///< for percolating cluster

    void PlotMatlab_DFN(string FileKey);
    //< matlab plot dfn

    void PlotMatlab_DFN_trim(string FileKey);
    //

    void PlotMatlab_DFN_and_Intersection(string FileKey);
    //

    void PlotMatlab_ORI_SCATTER(string FileKey);
    //

    void PlotMatlab_Traces_on_Model_surfaces(string FileKey);
    /////< Plot traces on surface

    void PlotMatlab_DFN_Highlight_Cluster(string FileKey);
    ///< Plot DFN highlighted by cluster values

    void PLotMatlab_DFN_Cluster_along_a_direction(string FileKey, string str);
    ///< Plot percolation cluster
    //spanning model along x, y or z axis

    void PlotMatlab_Radius_and_Area_kstest(string FileKey);
    //are both following the same distribution,
    // but seems not work

    void PlotMatlab_Radius_and_Perimeter(string FileKey);

    void DataFile_Radius_AreaAndPerimeter(string FileKey);
    ///< outputs the data

    void Create_whole_model_II(const Vector6d model_size, std::vector<std::vector<Vector3d>> Frac_verts);

    void Matlab_Out_Frac_matfile(string FileKey_mat);

    void Output_all_fractures_HDF5(string filename, string groupname, size_t field_width);
};

inline void Domain::Create_whole_model(const size_t n,
                                       const std::vector<double> DenWeight,
                                       gsl_rng *random_seed,
                                       const Vector6d model_size,
                                       const string str_ori,
                                       const string str_frac_size,
                                       const std::vector<Vector2d> array11,
                                       std::vector<Vector4d> array12,
                                       const std::vector<Vector7d> array13,
                                       const string conductivity_distri)
{

    Random_function r1 = Random_function(random_seed);

    Model_set(model_size);

    Last_frac_size = -1;

    if (mode_2D == false)
    {
        if (str_ori == "uniform")
        {
            for (size_t i = 0; i < n; ++i)
            {
                if (str_frac_size == "powerlaw")
                {
                    //cout << "debug 1\n";
                    Fracture f(str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);
                    //cout << "debug 2\n";
                    AddSquareFracture(i, f);

                    //cout << "debug 3\n";
                }
                else if (str_frac_size == "lognormal")
                {

                    Fracture f(str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);
                    AddSquareFracture(i, f);
                }
                else if (str_frac_size == "uniform")
                {
                    Fracture f(str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);
                    AddSquareFracture(i, f);
                }
                else if (str_frac_size == "single")
                {

                    Fracture f(str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);
                    AddSquareFracture(i, f);
                }
            }
        }
        else if (str_ori == "fisher")
        {
            size_t numofsets_1 = DenWeight.size();
            for (size_t i = 0; i < numofsets_1; ++i)
            {
                size_t init_jkk = 0;
                double Weight_k = 0;
                for (size_t j = 0; j < i; j++)
                {
                    Weight_k = Weight_k + DenWeight[j];
                }
                init_jkk = n * Weight_k;

                size_t end_jkk = 0;
                double Weight_i = 0;
                for (size_t j = 0; j <= i; ++j)
                {
                    Weight_i = DenWeight[j] + Weight_i;
                }
                end_jkk = n * Weight_i;
                //std::cout<<"start: "<< init_jkk<<" ; end: "<<end_jkk<<"\n";
                for (size_t j = init_jkk; j < end_jkk; ++j)
                {

                    if (str_frac_size == "powerlaw")
                    {
                        Fracture f(str_ori, str_frac_size, j, r1, array11, array12[i], array13[i], Last_frac_size, conductivity_distri);
                        AddSquareFracture(j, f);
                    }
                    else if (str_frac_size == "lognormal")
                    {
                        Fracture f(str_ori, str_frac_size, j, r1, array11, array12[i], array13[i], Last_frac_size, conductivity_distri);
                        AddSquareFracture(j, f);
                    }
                    else if (str_frac_size == "uniform")
                    {
                        Fracture f(str_ori, str_frac_size, j, r1, array11, array12[i], array13[i], Last_frac_size, conductivity_distri);
                        AddSquareFracture(j, f);
                    }
                    else if (str_frac_size == "single")
                    {
                        Fracture f(str_ori, str_frac_size, j, r1, array11, array12[i], array13[i], Last_frac_size, conductivity_distri);
                        AddSquareFracture(j, f);
                    }
                }
            }
        }
    }
    else
    {

        for (size_t i = 0; i < n; ++i)
        {
            bool mode2d = true;
            if (str_frac_size == "powerlaw")
            {
                //cout << "debug 1\n";
                Fracture f(mode2d, str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);
                //cout << "debug 2\n";
                AddSquareFracture(i, f);

                //cout << "debug 3\n";
            }
            else if (str_frac_size == "lognormal")
            {

                Fracture f(mode2d, str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);
                AddSquareFracture(i, f);
            }
            else if (str_frac_size == "uniform")
            {
                Fracture f(mode2d, str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);
                AddSquareFracture(i, f);
            }
            else if (str_frac_size == "single")
            {
                Fracture f(mode2d, str_ori, str_frac_size, i, r1, array11, array12[0] /*, array13*/, Last_frac_size, conductivity_distri);

                AddSquareFracture(i, f);
            }
        }
    }

    size_t nz = Fractures.size();
    if (nz == 0)
    {
        P32_total = 0;
        P32_connected = Eigen::VectorXd::Zero(3, 1);
        P30 = 0;
        P30_connected = Eigen::VectorXd::Zero(3, 1);

        Ratio_of_P32 = Eigen::VectorXd::Zero(3, 1);
        Ratio_of_P30 = Eigen::VectorXd::Zero(3, 1);

        return;
    };

    for (size_t i = 0; i < nz - 1; ++i)
    {
        for (size_t j = i + 1; j < nz; ++j)
        {
            if (mode_2D == false)
                Intersect(Fractures[i], Fractures[j], false);
            else
                Intersect(Fractures[i], Fractures[j], false);
        }
    }

    this->Clusters();

    this->Determine_max_cluster();
    this->Average_number_of_intersections_per_fracture();
};

inline void Domain::Model_set(const Vector6d model_size)
{
    double xmin = model_size[0];
    double xmax = model_size[1];
    double ymin = model_size[2];
    double ymax = model_size[3];
    double zmin = model_size[4];
    double zmax = model_size[5];
    /// top-----------
    size_t Tag_1 = 0;
    size_t Clus_1 = -1;

    std::vector<Vector3d> Verts_1;
    Verts_1.resize(4);
    Verts_1[0] << xmin, ymin, zmax;
    Verts_1[1] << xmax, ymin, zmax;
    Verts_1[2] << xmax, ymax, zmax;
    Verts_1[3] << xmin, ymax, zmax;

    Fracture Top(Tag_1, Clus_1, Verts_1);

    /// bottom-----------
    size_t Tag_2 = 1;
    size_t Clus_2 = -1;

    std::vector<Vector3d> Verts_2;
    Verts_2.resize(4);
    Verts_2[0] << xmin, ymin, zmin;
    Verts_2[1] << xmax, ymin, zmin;
    Verts_2[2] << xmax, ymax, zmin;
    Verts_2[3] << xmin, ymax, zmin;

    Fracture Bottom(Tag_2, Clus_2, Verts_2);

    /// front-----------
    size_t Tag_3 = 2;
    size_t Clus_3 = -1;

    std::vector<Vector3d> Verts_3;
    Verts_3.resize(4);
    Verts_3[0] << xmin, ymin, zmin;
    Verts_3[1] << xmax, ymin, zmin;
    Verts_3[2] << xmax, ymin, zmax;
    Verts_3[3] << xmin, ymin, zmax;

    Fracture Front(Tag_3, Clus_3, Verts_3);

    /// back-----------
    size_t Tag_4 = 3;
    size_t Clus_4 = -1;

    std::vector<Vector3d> Verts_4;
    Verts_4.resize(4);
    Verts_4[0] << xmin, ymax, zmin;
    Verts_4[1] << xmax, ymax, zmin;
    Verts_4[2] << xmax, ymax, zmax;
    Verts_4[3] << xmin, ymax, zmax;

    Fracture Back(Tag_4, Clus_4, Verts_4);

    /// left-----------
    size_t Tag_5 = 4;
    size_t Clus_5 = -1;

    std::vector<Vector3d> Verts_5;
    Verts_5.resize(4);
    Verts_5[0] << xmin, ymin, zmin;
    Verts_5[1] << xmin, ymax, zmin;
    Verts_5[2] << xmin, ymax, zmax;
    Verts_5[3] << xmin, ymin, zmax;

    Fracture Left(Tag_5, Clus_5, Verts_5);

    /// right-----------
    size_t Tag_6 = 5;
    size_t Clus_6 = -1;

    std::vector<Vector3d> Verts_6;
    Verts_6.resize(4);
    Verts_6[0] << xmax, ymin, zmin;
    Verts_6[1] << xmax, ymax, zmin;
    Verts_6[2] << xmax, ymax, zmax;
    Verts_6[3] << xmax, ymin, zmax;

    Fracture Right(Tag_6, Clus_6, Verts_6);

    ///--------push the six surfaces into Fractures, so now, remember, in Fractures, we really have (Size (of Fractures) minus six) fractures
    /// they are Fractures[0] to [5]
    Surfaces.push_back(Top);
    Surfaces.push_back(Bottom);
    Surfaces.push_back(Front);
    Surfaces.push_back(Back);
    Surfaces.push_back(Left);
    Surfaces.push_back(Right);

    Model_domain << zmax, zmin, ymin, ymax, xmin, xmax;
};

inline void Domain::AddSquareFracture(size_t Tag,
                                      Fracture &c)
{

    if (Model_domain(4) <= c.Center(0) &&
        Model_domain(5) >= c.Center(0) &&
        Model_domain(2) <= c.Center(1) &&
        Model_domain(3) >= c.Center(1) &&
        Model_domain(1) <= c.Center(2) &&
        Model_domain(0) >= c.Center(2))
    {
        bool y1 = Intersect_A(Surfaces[0], c);
        bool y2 = Intersect_A(Surfaces[1], c);
        bool y3 = Intersect_A(Surfaces[2], c);
        bool y4 = Intersect_A(Surfaces[3], c);
        bool y5 = Intersect_A(Surfaces[4], c);
        bool y6 = Intersect_A(Surfaces[5], c);
        if (y1 == 1 || y2 == 1 || y3 == 1 || y4 == 1 || y5 == 1 || y6 == 1)
        {

            if (y1 == 1)
            {
                c.If_intersect_surfaces(0) = 1;
            }

            if (y2 == 1)
            {
                c.If_intersect_surfaces(1) = 1;
            }
            if (y3 == 1)
            {
                c.If_intersect_surfaces(2) = 1;
            }
            if (y4 == 1)
            {
                c.If_intersect_surfaces(3) = 1;
            }
            if (y5 == 1)
            {
                c.If_intersect_surfaces(4) = 1;
            }
            if (y6 == 1)
            {
                c.If_intersect_surfaces(5) = 1;
            }

            std::vector<Vector3d> YT = c.Verts_trim;
            DFN::Intersection_between_polygon_and_3D_box Inse{YT, this->Model_domain};
            c.Verts_trim = YT;
            Modify_trimmed_fractures_attribute(c);
        }
        If_fracture_intersect_boundary(c);
        Fractures.push_back(c);
        Fractures[Fractures.size() - 1].Tag = Fractures.size() - 1;
        Last_frac_size = -1;
        c.Nvertices_trim = c.Verts_trim.size();
        return;
    }
    else
    {

        bool y1 = Intersect_A(Surfaces[0], c);
        bool y2 = Intersect_A(Surfaces[1], c);
        bool y3 = Intersect_A(Surfaces[2], c);
        bool y4 = Intersect_A(Surfaces[3], c);
        bool y5 = Intersect_A(Surfaces[4], c);
        bool y6 = Intersect_A(Surfaces[5], c);
        if (y1 == 1 || y2 == 1 || y3 == 1 || y4 == 1 || y5 == 1 || y6 == 1)
        {

            if (y1 == 1)
            {

                c.If_intersect_surfaces(0) = 1;
            }
            if (y2 == 1)
            {

                c.If_intersect_surfaces(1) = 1;
            }
            if (y3 == 1)
            {

                c.If_intersect_surfaces(2) = 1;
            }
            if (y4 == 1)
            {

                c.If_intersect_surfaces(3) = 1;
            }
            if (y5 == 1)
            {

                c.If_intersect_surfaces(4) = 1;
            }
            if (y6 == 1)
            {

                c.If_intersect_surfaces(5) = 1;
            }

            std::vector<Vector3d> YT = c.Verts_trim;
            DFN::Intersection_between_polygon_and_3D_box Inse{YT, this->Model_domain};

            if (YT.size() > 0)
            {
                c.Verts_trim = YT;
                If_fracture_intersect_boundary(c);
                Modify_trimmed_fractures_attribute(c);
                Fractures.push_back(c);
                Fractures[Fractures.size() - 1].Tag = Fractures.size() - 1;
                Last_frac_size = -1;
                c.Nvertices_trim = c.Verts_trim.size();
            }

            return;
        }
    }

    Last_frac_size = -1; //c.Radius;
    return;
};

inline bool Domain::Intersect_A(const Fracture F1, const Fracture F2)
{
    DFN::Polygon_convex_3D f1{F1.Verts};
    DFN::Polygon_convex_3D f2{F2.Verts};
    f1.Optimize();
    f2.Optimize();
    /*
    for (size_t i = 0; i < F1.Verts.size(); ++i)
        cout << F1.Verts[i].transpose() << endl;
    for (size_t i = 0; i < F2.Verts.size(); ++i)
        cout << F2.Verts[i].transpose() << endl;
    cout << endl;

    for (size_t i = 0; i < f1.Corners.size(); ++i)
        cout << f1.Corners[i].transpose() << endl;
    for (size_t i = 0; i < f2.Corners.size(); ++i)
        cout << f2.Corners[i].transpose() << endl;
    cout << endl;
    */
    //DFN::Intersection_Frac Interse{f1, f2};
    DFN::Intersection_Frac_boost Interse{f1, f2};
    return Interse.If_intersect;
}

inline void Domain::Modify_trimmed_fractures_attribute(Fracture &F2)
{
    F2.Nvertices_trim = F2.Verts_trim.size();
    ///Heron's formula
    F2.Area_trim = 0;
    for (size_t i = 0; i < F2.Nvertices_trim - 2; ++i)
    {
        size_t j = i + 1;

        size_t k = i + 2 - (size_t)((i + 2) / F2.Nvertices_trim) * (i + 2);

        double a, b, c, p;
        a = pow((F2.Verts_trim[0] - F2.Verts_trim[j]).dot((F2.Verts_trim[0] - F2.Verts_trim[j])), 0.5);
        b = pow((F2.Verts_trim[j] - F2.Verts_trim[k]).dot((F2.Verts_trim[j] - F2.Verts_trim[k])), 0.5);
        c = pow((F2.Verts_trim[k] - F2.Verts_trim[0]).dot((F2.Verts_trim[k] - F2.Verts_trim[0])), 0.5);
        p = (a + b + c) / 2;

        double Area_1;
        if (a == 0 || b == 0 || c == 0)
            Area_1 = 0;
        else
            Area_1 = pow((p * (p - a) * (p - b) * (p - c)), 0.5);

        F2.Area_trim = F2.Area_trim + Area_1;
    }

    F2.Perimeter_trim = 0;
    for (size_t i = 0; i < F2.Verts_trim.size(); ++i)
    {
        size_t j = i + 1 - (size_t)((i + 1) / (F2.Verts_trim.size())) * (i + 1);
        double p = pow((F2.Verts_trim[i] - F2.Verts_trim[j]).dot((F2.Verts_trim[i] - F2.Verts_trim[j])), 0.5);
        F2.Perimeter_trim = F2.Perimeter_trim + p;
    }
};

inline void Domain::If_fracture_intersect_boundary(Fracture &F2)
{
    if (F2.If_boundary.size() > 1)
    {
        throw Error_throw_pause("Error! In class of 'Domain', function 'If_fracture_intersect_boundary', fracture array 'If_boundary' should be initialized\n");
    }
    if (Surfaces.size() != 6)
    {
        throw Error_throw_pause("Error! In class of 'Domain', function 'If_fracture_intersect_boundary', Surfaces array 'Surfaces' should be initialized\n");
    }
    for (size_t i = 0; i < F2.Verts_trim.size(); ++i)
    {
        size_t ni = i + 1 - (size_t)((i + 1) / F2.Verts_trim.size()) * (i + 1);
        double top_surf = Surfaces[0].Verts[0](2);
        double bottom_surf = Surfaces[1].Verts[0](2);
        double front_surf = Surfaces[2].Verts[0](1);
        double back_surf = Surfaces[3].Verts[0](1);
        double left_surf = Surfaces[4].Verts[0](0);
        double right_surf = Surfaces[5].Verts[0](0);
        if (abs(F2.Verts_trim[i](2) - top_surf) < 0.001 && abs(F2.Verts_trim[ni](2) - top_surf))
        {
            Vector2d temoy;
            temoy << i, 0;
            F2.If_boundary.push_back(temoy);
        }
        if (abs(F2.Verts_trim[i](2) - bottom_surf) < 0.001 && abs(F2.Verts_trim[ni](2) - bottom_surf))
        {
            Vector2d temoy;
            temoy << i, 1;
            F2.If_boundary.push_back(temoy);
        }
        if (abs(F2.Verts_trim[i](1) - front_surf) < 0.001 && abs(F2.Verts_trim[ni](1) - front_surf))
        {
            Vector2d temoy;
            temoy << i, 2;
            F2.If_boundary.push_back(temoy);
        }
        if (abs(F2.Verts_trim[i](1) - back_surf) < 0.001 && abs(F2.Verts_trim[ni](1) - back_surf))
        {
            Vector2d temoy;
            temoy << i, 3;
            F2.If_boundary.push_back(temoy);
        }
        if (abs(F2.Verts_trim[i](0) - left_surf) < 0.001 && abs(F2.Verts_trim[ni](0) - left_surf))
        {
            Vector2d temoy;
            temoy << i, 4;
            F2.If_boundary.push_back(temoy);
        }
        if (abs(F2.Verts_trim[i](0) - right_surf) < 0.001 && abs(F2.Verts_trim[ni](0) - right_surf))
        {
            Vector2d temoy;
            temoy << i, 5;
            F2.If_boundary.push_back(temoy);
        }
    }
};

inline bool Domain::Intersect(Fracture &F1, Fracture &F2, const bool if_use_trim_frac)
{
    F1.Nvertices_trim = F1.Verts_trim.size();
    F2.Nvertices_trim = F2.Verts_trim.size();

    Vector3d A1;
    A1 << 0, 0, 0;

    Vector3d B1;
    B1 << 0, 0, 0;

    DFN::Polygon_convex_3D f1;
    DFN::Polygon_convex_3D f2;

    if (if_use_trim_frac == false)
    {
        DFN::Polygon_convex_3D x1{F1.Verts};
        DFN::Polygon_convex_3D x2{F2.Verts};
        f1.Create(x1);
        f2.Create(x2);
    }
    else
    {
        DFN::Polygon_convex_3D x1{F1.Verts_trim};
        DFN::Polygon_convex_3D x2{F2.Verts_trim};
        f1.Create(x1);
        f2.Create(x2);
    }

    f1.Optimize();
    f2.Optimize();

    //DFN::Intersection_Frac Interse{f1, f2};
    DFN::Intersection_Frac_boost Interse{f1, f2};
    //return Interse.If_intersect;
    if (Interse.If_intersect == true)
    {
        A1 << round(Interse.Intersection[0](0), 4),
            round(Interse.Intersection[0](1), 4),
            round(Interse.Intersection[0](2), 4);

        B1 << round(Interse.Intersection[1](0), 4),
            round(Interse.Intersection[1](1), 4),
            round(Interse.Intersection[1](2), 4);
    }
    else
    {
        return false;
    }

    ///std::cout<<"Intersection points are: \n"<<A1<<"\n"<<B1<<std::endl;

    //std::cout << "Found intersection between Fracture[" << F1.Tag << "] and Fracture[" << F2.Tag << "]!" << std::endl;
    ///Clusters(F1,F2, A1, B1);

    Connections.push_back(F1.Tag); //Fracture F1 and F2 are connected
    Connections.push_back(F2.Tag);
    std::pair<size_t, size_t> p = std::make_pair(F1.Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
    Vector3d x1, x2;
    x1 = A1;
    x2 = B1;
    Intersections[p] = std::make_pair(x1, x2);

    // std::set<size_t, size_t> Intersect_other_frac_after_trim;
    F1.Intersect_other_frac_after_trim.insert(F2.Tag);
    F2.Intersect_other_frac_after_trim.insert(F1.Tag);
    return true;
}

inline void Domain::Clusters()
{
    Graph GK(Fractures.size(), Connections);
    GK.CreateGraph_i(Listofclusters);

    //std::cout <<"Listofclusters.size(): " << Listofclusters.size() << std::endl;
    for (size_t i = 0; i < Listofclusters.size(); i++)
    {
        //std::cout << "debug 1.2" << std::endl;
        //std::cout<<Listofclusters[i]<<std::endl;
        for (size_t j = 0; j < Listofclusters[i].size(); j++)
        {
            //std::cout << "debug 1.3" << std::endl;
            Fractures[Listofclusters[i][j]].Clus = i;
        }
        //std::cout << "debug 1.4" << std::endl;
    }
}

inline void Domain::Determine_max_cluster()
{
    //Xi = 0;
    //max_R_s = 0;
    size_t max_cluster = 0;
    size_t max_size = 0;

    if (Listofclusters.size() != 0)
    {
        //max_cluster = 0;
        //max_size = 0;
        //double total_distance = 0;
        for (size_t i = 0; i < Listofclusters.size(); ++i)
        {
            //std::cout << "Cluster " << i << "'s size is: " << Listofclusters[i].size() << std::endl;
            if (max_size < Listofclusters[i].size())
            {
                max_size = Listofclusters[i].size();
                max_cluster = i;
            }
        }
        //std::cout << "the largest cluster is " << max_cluster << std::endl;
        //std::cout << "size of max cluster is " << Listofclusters[max_cluster].size() << std::endl;
        //cout << max_cluster << endl;
        //cout << Listofclusters[max_cluster].size() << endl;
        if (Listofclusters[max_cluster].size() > 1)
        {
            /*
            Center_of_cluster[0] = 0;
            Center_of_cluster[1] = 0;
            Center_of_cluster[2] = 0;

            for (size_t i = 0; i < Listofclusters[max_cluster].size(); ++i)
            {
                Center_of_cluster[0] = Center_of_cluster[0] + Fractures[Listofclusters[max_cluster][i]].Center[0];
                Center_of_cluster[1] = Center_of_cluster[1] + Fractures[Listofclusters[max_cluster][i]].Center[1];
                Center_of_cluster[2] = Center_of_cluster[2] + Fractures[Listofclusters[max_cluster][i]].Center[2];
            }
            Center_of_cluster[0] = Center_of_cluster[0] / Listofclusters[max_cluster].size();
            Center_of_cluster[1] = Center_of_cluster[1] / Listofclusters[max_cluster].size();
            Center_of_cluster[2] = Center_of_cluster[2] / Listofclusters[max_cluster].size();

            max_R_s = 0;
            for (size_t i = 0; i < Listofclusters[max_cluster].size(); ++i)
            {
                Vector3d distance = Fractures[Listofclusters[max_cluster][i]].Center - Center_of_cluster;
                double module = pow(pow(distance[0], 2) + pow(distance[1], 2) + pow(distance[2], 2), 0.5);
                if (max_R_s <= module)
                    max_R_s = module;
            }
            */
            double Model_volume = (Model_domain(0) - Model_domain(1)) * (Model_domain(3) - Model_domain(2)) * (Model_domain(5) - Model_domain(4));
            P30_largest_cluster = (double)Listofclusters[max_cluster].size() / Model_volume;
            //cout << "P30_largest_cluster: " << P30_largest_cluster << endl;
            double area_total = 0;
            for (size_t i = 0; i < Listofclusters[max_cluster].size(); ++i)
            {
                area_total = area_total + Fractures[Listofclusters[max_cluster][i]].Area;
            }

            P32_largest_cluster = area_total / Model_volume;
        }
        else
        {
            P30_largest_cluster = 0;
            P32_largest_cluster = 0;
        }
        //Xi = 0;
        //max_R_s = 0;
    }

    /// now, determine Xi
    /*
    if (Listofclusters.size() > 1)
    {
        size_t k = Listofclusters.size();
        size_t a[2][k];
        for (size_t i = 0; i < k; ++i)
        {
            a[0][i] = i;
            a[1][i] = 0;
        }

        std::vector<std::vector<size_t>> kss;
        kss.resize(max_size);

        for (size_t i = 0; i < max_size; ++i)
        {
            for (size_t j = 0; j < k; ++j)
            {
                if (a[1][j] == 0)
                {
                    if (Listofclusters[j].size() == i + 1)
                    {
                        kss[i].push_back(j);
                        a[1][j] = 1;
                    }
                }
            }
        }
        
        double numerator_h = 0;
        double denominator_h = 0;
        for (size_t i = 0; i < kss.size(); ++i)
        {
            for (size_t j = 0; j < kss[i].size(); ++j)
            {
                double R_s_squared = 0;
                if (i != 0)
                {
                    //std::cout << "\n----------------*\n";
                    double distance = 0;
                    for (size_t v = 0; v < Listofclusters[kss[i][j]].size(); ++v)
                    {
                        for (size_t w = 0; w < Listofclusters[kss[i][j]].size(); ++w)
                        {

                            double auu = Fractures[Listofclusters[kss[i][j]][v]].Center.dot(Fractures[Listofclusters[kss[i][j]][v]].Center) + Fractures[Listofclusters[kss[i][j]][w]].Center.dot(Fractures[Listofclusters[kss[i][j]][w]].Center) - 2 * Fractures[Listofclusters[kss[i][j]][v]].Center.dot(Fractures[Listofclusters[kss[i][j]][w]].Center);
                            distance = distance + auu;

                            //std::cout << "\tdistance between Fracture " << Listofclusters[kss[i][j]][v] << " and Fracture " << Listofclusters[kss[i][j]][w] << " is " << auu << ";\n";
                        }
                    }
                    //R_s_squared = 0.5 * distance / (double)pow((i + 1), (i + 1));
                    R_s_squared = 0.5 * distance / (double)pow((i + 1), 2);
                    //std::cout << "Counting distance of cluster " << kss[i][j] << ", and it's size is (" << i + 1 << "): ";
                    //std::cout << R_s_squared << "\n----------------\n";
                }
                double numerator_s = R_s_squared * (i + 1) * (i + 1) * ((double)kss[i].size() / Fractures.size());
                double denominator_s = (i + 1) * (i + 1) * ((double)kss[i].size() / Fractures.size());

                numerator_h = numerator_h + numerator_s;
                denominator_h = denominator_h + denominator_s;
            }
        }
        //std::cout << numerator_h <<"\n" << denominator_h <<"\n";
        Xi = 2 * numerator_h / denominator_h;
        Xi = pow(Xi, 0.5);*/
    /*
        Xi = 0;
        size_t ku = 0;
        for (size_t i = 1; i < kss.size(); ++i)
        {
            if (kss[i].size() != 0)
            {
                ku++;
                double D_a = 0;
                for (size_t j = 0; j < kss[i].size(); ++j)
                {
                    double D_c = 0;
                    for (size_t v = 0; v < Listofclusters[kss[i][j]].size(); ++v)
                    {
                        size_t ID_frac = Listofclusters[kss[i][j]][v];
                        //std::vector<size_t> Connections
                        std::vector<size_t> connected_frac;
                        for (size_t w = 0; w < Connections.size(); w++)
                        {
                            if (Connections[w] == ID_frac)
                            {
                                if (w % 2 == 0)
                                    connected_frac.push_back(Connections[w + 1]);
                                else if (w % 2 == 1)
                                    connected_frac.push_back(Connections[w - 1]);
                            }
                        }

                        double dist = 0;
                        for (size_t w = 0; w < connected_frac.size(); ++w)
                        {
                            dist += (Fractures[ID_frac].Center - Fractures[connected_frac[w]].Center).norm();
                        }

                        dist = dist / connected_frac.size();
                        //cout << "dist: " << dist << endl;
                        D_c += dist;
                    }
                    D_c = D_c / Listofclusters[kss[i][j]].size();
                    D_a += D_c;
                    //cout << "D_c: " << D_c << endl;
                }
                D_a = D_a / kss[i].size();
                //cout << "D_a: " << D_a << endl;
                Xi += D_a;
            }
        }
        if (ku != 0)
            Xi = Xi / ku;
        else
            Xi = 0;
    }
    else
    {
        Xi = 0;
    }
    */
    //std::cout << Center_of_cluster[0] << ", " << Center_of_cluster[1] <<", " << Center_of_cluster[2] << std::endl;
};

inline void Domain::Average_number_of_intersections_per_fracture()
{
    n_I = (double)(Connections.size() / 2) / (double)Fractures.size();
};

inline void Domain::Determine_excluded_volume(const string str_ori, const string str_frac_size, double alpha_g, double kappa, double mean_i, double var_i, double min_R_i, double max_R_i){
    //
    //if (str_ori == "uniform")
    //{
    //    //double C1 = (1-alpha_g)/((2-alpha_g)*(pow(x1,1-alpha_g)-pow(x0,1-alpha_g)));
    //    //double Expected_value_of_radius = C1*pow(x1,2-alpha_g) - C1*pow(x0,2-alpha_g);
    //    //std::cout<<"Expected_value_of_radius: "<<Expected_value_of_radius<<std::endl;
    //    if (str_frac_size == "powerlaw")
    //    {
    //        double x0, x1;
    //        x0 = min_R_i * pow(2.0, 0.5);
    //        x1 = max_R_i * pow(2.0, 0.5);
    //        double C2 = (1 - alpha_g) / (pow(x1, 1 - alpha_g) - pow(x0, 1 - alpha_g));
    //        //********below is 0.5*<A*P>
    //        Excluded_volume_1 = C2 * 4.0 * (pow(x1, 4.0 - alpha_g) - pow(x0, 4.0 - alpha_g)) / (4.0 - alpha_g);
    //        Excluded_volume_1 = 0.5 * Excluded_volume_1;
    //
    //        //********below is 0.5*<A>*<P>
    //        Excluded_volume_2 = 0.5 * C2 * (pow(x1, 3.0 - alpha_g) - pow(x0, 3.0 - alpha_g)) / (3.0 - alpha_g);                       //0.5*<A>
    //        Excluded_volume_2 = Excluded_volume_2 * (4.0 * C2 * (pow(x1, 2.0 - alpha_g) - pow(x0, 2.0 - alpha_g)) / (2.0 - alpha_g)); //<A>*<P>
    //
    //        //********below is 0.5 * <l^3>
    //        Excluded_volume_3 = 0.5 * C2 * (pow(x1, 4.0 - alpha_g) - pow(x0, 4.0 - alpha_g)) / (4.0 - alpha_g);
    //
    //        //********below is 0.5 * <l^3>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (C2 / (3.0 - alpha_g) * (pow(x1, 3.0 - alpha_g) - pow(x0, 3.0 - alpha_g)));
    //
    //        //********below is 0.5 * <A*P>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (C2 / (3.0 - alpha_g) * (pow(x1, 3.0 - alpha_g) - pow(x0, 3.0 - alpha_g)));
    //    }
    //    else if (str_frac_size == "lognormal")
    //    {
    //        double mean, var;
    //        mean = mean_i * pow(2, 0.5);
    //        var = 2.0 * var_i;
    //
    //        double mean_1 = log(mean * mean / (pow(var + mean * mean, 0.5)));
    //        double var_1 = pow(log(1 + ((double)var) / (mean * mean)), 0.5); //var_1 is input std. deviation
    //
    //        //********below is 0.5*<A*P>
    //        Excluded_volume_1 = 0.5 * 4.0 * exp(3 * mean_1 + (9.00 / 2.0) * (var_1 * var_1));
    //        //std::cout << "Vex: " << Excluded_volume_1 <<"\n";
    //
    //        //********below is 0.5*<A>*<P>
    //        Excluded_volume_2 = 0.5 * exp(2.0 * mean_1 + 2.0 * var_1 * var_1) * 4.0 * exp(mean_1 + 0.5 * var_1 * var_1);
    //
    //        //********below is 0.5 * <l^3>
    //        Excluded_volume_3 = 0.5 * exp(3.0 * mean_1 + (9.00 / 2.0) * (var_1 * var_1));
    //
    //        //********below is 0.5 * <l^3>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (exp(2.0 * mean_1 + 2.0 * var_1 * var_1));
    //
    //        //********below is 0.5 * <A*P>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (exp(2.0 * mean_1 + 2.0 * var_1 * var_1));
    //    }
    //    else if (str_frac_size == "uniform")
    //    {
    //        double max_R, min_R;
    //        max_R = max_R_i * pow(2.0, 0.5);
    //        min_R = min_R_i * pow(2.0, 0.5);
    //        double C2 = 1 / (max_R - min_R);
    //
    //        //********below is 0.5*<A*P>
    //        Excluded_volume_1 = 0.5 * C2 * 4.0 / 4.0 * (pow(max_R, 4) - pow(min_R, 4));
    //
    //        //********below is 0.5*<A>*<P>
    //        Excluded_volume_2 = 0.5 * C2 * (pow(max_R, 3.0) - pow(min_R, 3.0)) / 3.0 * C2 * 4.0 * (pow(max_R, 2.0) - pow(min_R, 2.0)) / 2.0;
    //
    //        //********below is 0.5 * <l^3>
    //        Excluded_volume_3 = 0.5 * C2 * ((pow(max_R, 4.0) - pow(min_R, 4.0)) / 4.0);
    //
    //        //********below is 0.5 * <l^3>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (C2 * ((pow(max_R, 3.0) - pow(min_R, 3.0)) / 3.0));
    //
    //        //********below is 0.5 * <A*P>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (C2 * ((pow(max_R, 3.0) - pow(min_R, 3.0)) / 3.0));
    //    }
    //    else if (str_frac_size == "single")
    //    {
    //        double R = min_R_i * pow(2.0, 0.5);
    //        //********below is 0.5*<A*P>
    //        Excluded_volume_1 = 0.5 * pow(R, 2) * 4.0 * R;
    //
    //        //********below is 0.5*<A>*<P>
    //        Excluded_volume_2 = Excluded_volume_1;
    //
    //        //********below is 0.5 * <l^3>
    //        Excluded_volume_3 = 0.5 * R * R * R;
    //
    //        //********below is 0.5 * <l^3>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (R * R);
    //
    //        //********below is 0.5 * <A*P>/<l^2>, and note that the percolation parameter equals to P32*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (R * R);
    //    }
    //}
    //else if (str_ori == "fisher")
    //{
    //    //here, we need to calculate more,i.e., the statistical average of sinγ
    //    double psi = (2 / ((sinh(kappa)) * (sinh(kappa)))) * (first_modified_Bessel(0, 2 * kappa) - (1 / kappa) * first_modified_Bessel(1, 2 * kappa));
    //    double pi = acos(-1);
    //    double SinGamma = psi * pi / 4;
    //    //--------------------------------------------------------
    //    if (str_frac_size == "powerlaw")
    //    {
    //        double x0, x1;
    //        x0 = min_R_i * pow(2.0, 0.5);
    //        x1 = max_R_i * pow(2.0, 0.5);
    //        double C2 = (1 - alpha_g) / (pow(x1, 1 - alpha_g) - pow(x0, 1 - alpha_g));
    //        //********below is 2/pi * < sinγ > * < A * P >
    //        Excluded_volume_1 = C2 * 4.0 * (pow(x1, 4.0 - alpha_g) - pow(x0, 4.0 - alpha_g)) / (4.0 - alpha_g);
    //        Excluded_volume_1 = 2 / pi * SinGamma * Excluded_volume_1;
    //
    //        //********below is 2/pi * < sinγ > *<A>*<P>
    //        Excluded_volume_2 = C2 * (pow(x1, 3.0 - alpha_g) - pow(x0, 3.0 - alpha_g)) / (3.0 - alpha_g);                                                 //<A>
    //        Excluded_volume_2 = 2 / pi * SinGamma * Excluded_volume_2 * (4.0 * C2 * (pow(x1, 2.0 - alpha_g) - pow(x0, 2.0 - alpha_g)) / (2.0 - alpha_g)); //<A>*<P>
    //
    //        //********below is 2/pi * < sinγ > * <l^3>
    //        Excluded_volume_3 = 2 / pi * SinGamma * C2 * (pow(x1, 4.0 - alpha_g) - pow(x0, 4.0 - alpha_g)) / (4.0 - alpha_g);
    //
    //        //********below is 2/pi * < sinγ > * <l^3>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (C2 / (3.0 - alpha_g) * (pow(x1, 3.0 - alpha_g) - pow(x0, 3.0 - alpha_g)));
    //
    //        //********below is 2/pi * < sinγ > * <A*P>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (C2 / (3.0 - alpha_g) * (pow(x1, 3.0 - alpha_g) - pow(x0, 3.0 - alpha_g)));
    //    }
    //    else if (str_frac_size == "lognormal")
    //    {
    //        double mean, var;
    //        mean = mean_i * pow(2, 0.5);
    //        var = 2.0 * var_i;
    //
    //        double mean_1 = log(mean * mean / (pow(var + mean * mean, 0.5)));
    //        double var_1 = pow(log(1 + ((double)var) / (mean * mean)), 0.5); //var_1 is input std. deviation
    //
    //        //********below is 2/pi * < sinγ > * < A * P >
    //        Excluded_volume_1 = 2 / pi * SinGamma * 4.0 * exp(3 * mean_1 + (9.00 / 2.0) * (var_1 * var_1));
    //        //std::cout << "Vex: " << Excluded_volume_1 <<"\n";
    //
    //        //********below is 2/pi * < sinγ > *<A>*<P>
    //        Excluded_volume_2 = 2 / pi * SinGamma * exp(2.0 * mean_1 + 2.0 * var_1 * var_1) * 4.0 * exp(mean_1 + 0.5 * var_1 * var_1);
    //
    //        //********below is 2/pi * < sinγ > * <l^3>
    //        Excluded_volume_3 = 2 / pi * SinGamma * exp(3.0 * mean_1 + (9.00 / 2.0) * (var_1 * var_1));
    //
    //        //********below is 2/pi * < sinγ > * <l^3>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (exp(2.0 * mean_1 + 2.0 * var_1 * var_1));
    //
    //        //********below is 2/pi * < sinγ > * <A*P>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (exp(2.0 * mean_1 + 2.0 * var_1 * var_1));
    //    }
    //    else if (str_frac_size == "uniform")
    //    {
    //        double max_R, min_R;
    //        max_R = max_R_i * pow(2.0, 0.5);
    //        min_R = min_R_i * pow(2.0, 0.5);
    //        double C2 = 1 / (max_R - min_R);
    //
    //        //********below is 2/pi * < sinγ > * < A * P >
    //        Excluded_volume_1 = 2 / pi * SinGamma * C2 * 4.0 / 4.0 * (pow(max_R, 4) - pow(min_R, 4));
    //
    //        //********below is 2/pi * < sinγ > *<A>*<P>
    //        Excluded_volume_2 = 2 / pi * SinGamma * C2 * (pow(max_R, 3.0) - pow(min_R, 3.0)) / 3.0 * C2 * 4.0 * (pow(max_R, 2.0) - pow(min_R, 2.0)) / 2.0;
    //
    //        //********below is 2/pi * < sinγ > * <l^3>
    //        Excluded_volume_3 = 2 / pi * SinGamma * C2 * ((pow(max_R, 4.0) - pow(min_R, 4.0)) / 4.0);
    //
    //        //********below is 2/pi * < sinγ > * <l^3>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (C2 * ((pow(max_R, 3.0) - pow(min_R, 3.0)) / 3.0));
    //
    //        //********below is 2/pi * < sinγ > * <A*P>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (C2 * ((pow(max_R, 3.0) - pow(min_R, 3.0)) / 3.0));
    //    }
    //    else if (str_frac_size == "single")
    //    {
    //        double R = min_R_i * pow(2.0, 0.5);
    //        //********below is 2/pi * < sinγ > * < A * P >
    //        Excluded_volume_1 = 2 / pi * SinGamma * (R * R) * (4.0 * R);
    //
    //        //********below is 2/pi * < sinγ > *<A>*<P>
    //        Excluded_volume_2 = Excluded_volume_1;
    //
    //        //********below is 2/pi * < sinγ > * <l^3>
    //        Excluded_volume_3 = 2 / pi * SinGamma * R * R * R;
    //
    //        //********below is 2/pi * < sinγ > * <l^3>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_4 = Excluded_volume_3 / (R * R);
    //
    //        //********below is 2/pi * < sinγ > * <A*P>/<l^2>, and note that the percolation parameter equals to P32*< sinγ >*<l^3>/<l^2>
    //        Excluded_volume_5 = Excluded_volume_1 / (R * R);
    //    }
    //}
    //else
    //{
    //    throw Error_throw_pause("Error! Please define the orientation distribution!\n");
    //}
    //
};

inline void Domain::Identify_percolation_clusters()
{
    Percolation_cluster.resize(3);

    if (Fractures.size() == 0)
    {
        //std::cout << "debug!!!?\n";
        return;
    }
    for (size_t i = 0; i < Listofclusters.size(); ++i) //Z direction
    {
        size_t q1 = 0, q2 = 0;
        for (size_t j = 0; j < Listofclusters[i].size(); ++j)
        {
            if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(0) == 1)
                q1 = 1;
            if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(1) == 1)
                q2 = 1;
            if (q1 == 1 && q2 == 1)
            {
                Percolation_cluster[2].push_back(i);
                break;
            }
        }
    }

    for (size_t i = 0; i < Listofclusters.size(); ++i) //Y direction
    {
        size_t q1 = 0, q2 = 0;
        for (size_t j = 0; j < Listofclusters[i].size(); ++j)
        {
            if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(2) == 1)
                q1 = 1;
            if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(3) == 1)
                q2 = 1;
            if (q1 == 1 && q2 == 1)
            {
                Percolation_cluster[1].push_back(i);
                break;
            }
        }
    }

    for (size_t i = 0; i < Listofclusters.size(); ++i) //X direction
    {
        size_t q1 = 0, q2 = 0;
        for (size_t j = 0; j < Listofclusters[i].size(); ++j)
        {
            if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(4) == 1)
                q1 = 1;
            if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(5) == 1)
                q2 = 1;
            if (q1 == 1 && q2 == 1)
            {
                Percolation_cluster[0].push_back(i);
                break;
            }
        }
    }

    if (Percolation_cluster[0].size() > 0)
    {
        Percolation_status[0] = true;
    }

    if (Percolation_cluster[1].size() > 0)
    {
        Percolation_status[1] = true;
    }

    if (Percolation_cluster[2].size() > 0)
    {
        Percolation_status[2] = true;
    }
}

inline void Domain::Connectivity_analysis()
{

    for (size_t km = 0; km < Percolation_cluster.size(); ++km)
    {
        std::vector<size_t> temp;
        temp = Percolation_cluster[km];

        P32_connected[km] = 0;
        P32_total = 0;
        double Area_connected = 0;
        double Area_total = 0;
        double Model_volume = (Model_domain(0) - Model_domain(1)) * (Model_domain(3) - Model_domain(2)) * (Model_domain(5) - Model_domain(4));

        double nooffractures_connected = 0;

        // temp is the percolating cluster
        for (size_t i = 0; i < temp.size(); ++i)
        {
            for (size_t j = 0; j < Listofclusters[temp[i]].size(); ++j)
            {
                size_t nf = Listofclusters[temp[i]][j];
                nooffractures_connected = nooffractures_connected + 1;
                Area_connected = Area_connected + Fractures[nf].Area;
            }
        }
        P32_connected[km] = Area_connected / Model_volume;
        for (size_t i = 0; i < Fractures.size(); ++i)
        {
            Area_total = Area_total + Fractures[i].Area;
        }

        P32_total = Area_total / Model_volume;
        Ratio_of_P32[km] = P32_connected[km] / P32_total;

        P30 = Fractures.size() / Model_volume;
        P30_connected[km] = nooffractures_connected / Model_volume;

        Ratio_of_P30[km] = P30_connected[km] / P30;
    }
}; // namespace DFN

inline void Domain::PlotMatlab_DFN(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);

    //Plotting the fractures
    oss << "clc;\nclose all;\nclear all;\n";
    for (size_t nf = 0; nf < Fractures.size(); nf++)
    {
        size_t n_verts = Fractures[nf].Verts.size();
        oss << "frac" << nf + 1 << " = fill3([";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts[nv_1](0) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts[nv_1](1) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts[nv_1](2) << " ";
        }
        oss << "],[rand rand rand]);\ngrid on;\nhold on;\n";
    }

    //Plotting the model domain
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
        }
    }
    double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
    double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
    double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

    oss.close();
};

inline void Domain::PlotMatlab_DFN_trim(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    //Plotting the fractures
    for (size_t nf = 0; nf < Fractures.size(); nf++)
    {
        size_t n_verts = Fractures[nf].Verts_trim.size();
        oss << "fill3([";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts_trim[nv_1](0) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts_trim[nv_1](1) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts_trim[nv_1](2) << " ";
        }
        oss << "],[rand rand rand]);\ngrid on;\nhold on;\n";
    }

    //Plotting the model domain
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
        }
    }
    double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
    double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
    double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

    //Open Matlab script to plot
    oss.close();
};

inline void Domain::PlotMatlab_DFN_and_Intersection(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);

    oss << "clc;\nclose all;\nclear all;\n";
    //Plotting the fractures
    for (size_t nf = 0; nf < Fractures.size(); nf++)
    {
        size_t n_verts = Fractures[nf].Verts.size();
        oss << "frac" << nf + 1 << " = fill3([";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts[nv_1](0) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts[nv_1](1) << " ";
        }
        oss << "],[";
        for (size_t nv = 0; nv < n_verts + 1; ++nv)
        {
            size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
            oss << Fractures[nf].Verts[nv_1](2) << " ";
        }
        oss << "],[rand rand rand]);\ngrid on;\nhold on;\n";
    }

    //Plotting the intersections
    //std::cout<<Connections.size()<<std::endl;

    for (std::map<std::pair<size_t, size_t>, std::pair<Vector3d, Vector3d>>::iterator its = Intersections.begin();
         its != Intersections.end(); its++)
    {
        size_t i1 = its->first.first;
        size_t i2 = its->first.second;
        Vector3d x1 = Intersections[std::make_pair(i1, i2)].first;
        Vector3d x2 = Intersections[std::make_pair(i1, i2)].second;
        if ((x1 - x2).norm() > 0.01)
        {
            oss << "plot3(";
            oss << "[" << x1(0) << " " << x2(0) << "],";
            oss << "[" << x1(1) << " " << x2(1) << "],";
            oss << "[" << x1(2) << " " << x2(2) << "],";
            oss << "'color',[0 0 1],'Linewidth',3);\ngrid on;\nhold on;\n";
        }
        else
        {
            oss << "scatter3(" << x1(0) << "," << x1(1) << ", " << x1(2) << ", '*', 'linewidth', 3);\nhold on;\n";
        }
    }

    //Plotting the model domain
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
        }
    }
    double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
    double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
    double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

    oss.close();
}

inline void Domain::PlotMatlab_ORI_SCATTER(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    double pi = acos(-1);
    oss << "th = [";
    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        double DD = Fractures[i].Dip_direction;

        double alpha = 0;
        if (DD > 90)
            alpha = 450 - DD;
        else if (DD <= 90)
            alpha = 90 - DD;
        alpha = alpha * pi / 180.0;
        oss << alpha << ", ...\n";

        //oss << DD * pi / 180.0 << ", ...\n";
    }
    oss << "];\nr = [";

    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        double DA = Fractures[i].Dip_angle;
        double beta = DA;
        beta = beta * pi / 180.0;
        oss << beta << ", ...\n";
    }
    oss << "];\npolarscatter(th,r,'filled');\nhold on;\nrlim([0 " << pi / 2 << "]);\n";
    oss << "hold on;\nrticks([" << pi / 12 << " " << 2 * pi / 12 << " " << 3 * pi / 12 << " " << 4 * pi / 12 << " " << 5 * pi / 12 << " " << 6 * pi / 12 << "]);\n";
    oss << "set(gca,'thetaticklabel',[]);\n";
    oss << "set(gca,'rticklabel',[]);";
    oss.close();
};

inline void Domain::PlotMatlab_Traces_on_Model_surfaces(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    //plotting the intersections between model surfaces and fractures

    for (size_t nc = 0; nc < Connections_S.size() / 2; nc++)
    {
        size_t i1 = Connections_S[2 * nc];
        size_t i2 = Connections_S[2 * nc + 1];
        Vector3d x1 = Intersections_S[std::make_pair(i1, i2)].first;
        Vector3d x2 = Intersections_S[std::make_pair(i1, i2)].second;
        oss << "plot3(";
        oss << "[" << x1(0) << " " << x2(0) << "],";
        oss << "[" << x1(1) << " " << x2(1) << "],";
        oss << "[" << x1(2) << " " << x2(2) << "],";
        oss << "'color',[0 0 1],'Linewidth',3);\ngrid on;\nhold on;\n";
    }

    //Plotting the model domain
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
        }
    }
    double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
    double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
    double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";
    oss.close();
};

inline void Domain::PlotMatlab_DFN_Highlight_Cluster(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    //Plotting the fractures, they are distinguished by clus values
    for (size_t i = 0; i < Listofclusters.size(); ++i)
    {
        double rand_1 = random_double(0, 1);
        double rand_2 = random_double(0, 1);
        double rand_3 = random_double(0, 1);
        for (size_t j = 0; j < Listofclusters[i].size(); ++j)
        {
            size_t nf = Listofclusters[i][j];
            size_t n_verts = Fractures[nf].Verts.size();
            oss << "fill3([";
            for (size_t nv = 0; nv < n_verts + 1; ++nv)
            {
                size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                oss << Fractures[nf].Verts[nv_1](0) << " ";
            }
            oss << "],[";
            for (size_t nv = 0; nv < n_verts + 1; ++nv)
            {
                size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                oss << Fractures[nf].Verts[nv_1](1) << " ";
            }
            oss << "],[";
            for (size_t nv = 0; nv < n_verts + 1; ++nv)
            {
                size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                oss << Fractures[nf].Verts[nv_1](2) << " ";
            }
            oss << "],[" << rand_1 << " " << rand_2 << " " << rand_3 << "]);\ngrid on;\nhold on;\n";
        }
    }

    //Plotting the model domain
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
        }
    }
    double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
    double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
    double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

    oss.close();
};

inline void Domain::PLotMatlab_DFN_Cluster_along_a_direction(string FileKey, string str)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    //Plotting the fractures
    std::vector<size_t> temp;
    if (str == "x")
    {
        temp.resize(Percolation_cluster[0].size());
        temp = Percolation_cluster[0];
    }
    else if (str == "y")
    {
        temp.resize(Percolation_cluster[1].size());
        temp = Percolation_cluster[1];
    }
    else if (str == "z")
    {
        temp.resize(Percolation_cluster[2].size());
        temp = Percolation_cluster[2];
    }
    else
    {
        throw Error_throw_pause("please define the percolation direction with char 'x', 'y' or 'z'\n");
    };

    for (size_t i = 0; i < temp.size(); ++i)
    {
        double rand_1 = random_double(0, 1);
        double rand_2 = random_double(0, 1);
        double rand_3 = random_double(0, 1);
        for (size_t j = 0; j < Listofclusters[temp[i]].size(); ++j)
        {
            size_t nf = Listofclusters[temp[i]][j];
            size_t n_verts = Fractures[nf].Verts.size();
            oss << "fill3([";
            for (size_t nv = 0; nv < n_verts + 1; ++nv)
            {
                size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                oss << Fractures[nf].Verts[nv_1](0) << " ";
            }
            oss << "],[";
            for (size_t nv = 0; nv < n_verts + 1; ++nv)
            {
                size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                oss << Fractures[nf].Verts[nv_1](1) << " ";
            }
            oss << "],[";
            for (size_t nv = 0; nv < n_verts + 1; ++nv)
            {
                size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
                oss << Fractures[nf].Verts[nv_1](2) << " ";
            }
            oss << "],[" << rand_1 << " " << rand_2 << " " << rand_3 << "]);\ngrid on;\nhold on;\n";
        }
    }

    //Plotting the model domain
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
            oss << "plot3(";
            oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
            oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
            oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
            oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
        }
    }
    double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
    double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
    double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
    oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";
    oss.close();
};

inline void Domain::PlotMatlab_Radius_and_Area_kstest(string FileKey)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "x1=[";
    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        oss << Fractures[i].Radius << "\t";
    }
    oss << "];\n";

    oss << "x2=[";
    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        oss << Fractures[i].Area << "\t";
    }
    oss << "];\n";

    oss << "[h,p,k] = kstest2(x1,x2);\n";
    oss << "%if h = 1, means the two groups of data are not having similar distributions;\n";
    oss << "hold on;\n";
    oss << "nbins = 30;\n";
    oss << "subplot(2,1,1);\n";
    oss << "histogram(x1,nbins);\n";
    oss << "hold on;\n";
    oss << "subplot(2,1,2);\n";
    oss << "%histogram(x2);\n";
    oss << "histogram(x2,nbins);\n";

    oss.close();
};

inline void Domain::PlotMatlab_Radius_and_Perimeter(string FileKey)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "Radius=[";
    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        oss << Fractures[i].Radius << ", ...\n";
    }
    oss << "];\n";
    /*
    oss << "x2=[";
    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        oss << Fractures[i].Perimeter << "\t";
    }
    oss << "];\n";
    */
    oss << "nbins = 30;\n";
    oss << "figure(1);\n";
    oss << "histogram(Radius,nbins);\n";
    oss << "hold on;\n";
    oss << "[Frequency, Data_bin]=hist(Radius,nbins);\n";
    oss << "\n\n%if power law, uncommect the following-------\n";
    oss << "%f = fittype('c*x^(-alpha)','independent','x','coefficients',{'c','alpha'} );\n";
    oss << "%cfun = fit(Data_bin',Frequency',f)\n";
    oss << "%upper = ;\n";
    oss << "%lower = ;\n";
    oss << "%alph = ;\n";
    oss << "%true_constant = (1-alph)/(upper^(1-alph) - lower^(1-alph))\n\n\n";

    oss << "\n\n%if lognormal, uncommect the following-------\n";
    oss << "%figure(2);\n";
    oss << "%histfit(Radius, nbins, 'lognormal');\n";
    oss << "%pd=fitdist(Radius','lognormal')\n";
    oss << "%Ex = exp(pd.mu + pd.sigma^2*0.5)\n";
    oss << "%Dx = exp(2*pd.mu+pd.sigma^2)*(exp(pd.sigma^2)-1)\n";
    oss.close();
};

inline void Domain::DataFile_Radius_AreaAndPerimeter(string FileKey)
{
    std::ofstream oss(FileKey, ios::out);
    oss << "Radius = [";
    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        oss << Fractures[i].Radius << ";\n";
    };
    oss << "];\n";
    oss.close();
}

inline void Domain::Re_identify_intersection_considering_trimmed_frac()
{
    Connections.clear();
    Listofclusters.clear();
    Percolation_cluster.clear();
    Intersections.erase(Intersections.begin(), Intersections.end());
    Percolation_status = {false, false, false};

    size_t nz = Fractures.size();

    for (size_t i = 0; i < nz - 1; ++i)
    {
        for (size_t j = i + 1; j < nz; ++j)
        {
            Intersect(Fractures[i], Fractures[j], true);
        }
    }

    //#pragma omp critical
    //{
    Clusters();
    //}
};

inline void Domain::Create_whole_model_II(const Vector6d model_size, std::vector<std::vector<Vector3d>> Frac_verts)
{
    this->Model_set(model_size);

    for (size_t i = 0; i < Frac_verts.size(); ++i)
    {
        Fracture f(i, -10, Frac_verts[i]);

        AddSquareFracture(i, f);
    }

    size_t nz = Fractures.size();

    for (size_t i = 0; i < nz - 1; ++i)
    {
        for (size_t j = i + 1; j < nz; ++j)
        {
            Intersect(Fractures[i], Fractures[j], false);
        }
    }

    //#pragma omp critical
    //  {
    Clusters();
    //}
}

inline void Domain::Matlab_Out_Frac_matfile(string FileKey_mat)
{
    for (size_t i = 0; i < this->Fractures.size(); ++i)
    {
        size_t len = Fractures[i].Verts.size();

        vector<double> pData1(len), pData2(len), pData3(len);

        for (size_t j = 0; j < len; j++)
        {
            pData1[j] = Fractures[i].Verts[j](0);
            pData2[j] = Fractures[i].Verts[j](1);
            pData3[j] = Fractures[i].Verts[j](2);
        }

        string ft = to_string(i + 1);

        string Fracx = "Frac_" + ft + "_x";
        string Fracy = "Frac_" + ft + "_y";
        string Fracz = "Frac_" + ft + "_z";

        DFN::MATLAB_DATA_API M_1;

        if (i == 0)
            M_1.Write_mat(FileKey_mat, "w", len, len, 1, pData1, Fracx);
        else
            M_1.Write_mat(FileKey_mat, "u", len, len, 1, pData1, Fracx);
        M_1.Write_mat(FileKey_mat, "u", len, len, 1, pData2, Fracy);
        M_1.Write_mat(FileKey_mat, "u", len, len, 1, pData3, Fracz);
    }

    DFN::MATLAB_DATA_API M_1;
    M_1.Write_mat(FileKey_mat, "u", 1, 1, 1, vector<double>{(double)this->Fractures.size()}, "Num_fracs");
};

inline void Domain::Output_all_fractures_HDF5(string filename, string groupname, size_t field_width)
{
    DFN::HDF5_DFN h5file;
    vector<string> string_field_overall;
    vector<vector<double>> data_set_overall;

    for (size_t i = 0; i < Fractures.size(); ++i)
    {
        string Frac_name = "Frac_" + To_string_with_width(i + 1, field_width) + "_";

        vector<double> Verts_x(Fractures[i].Verts.size()),
            Verts_y(Fractures[i].Verts.size()),
            Verts_z(Fractures[i].Verts.size());
        for (size_t j = 0; j < Fractures[i].Verts.size(); ++j)
        {
            Verts_x[j] = Fractures[i].Verts[j][0];
            Verts_y[j] = Fractures[i].Verts[j][1];
            Verts_z[j] = Fractures[i].Verts[j][2];
        }

        vector<double> Verts_trim_x(Fractures[i].Verts_trim.size()),
            Verts_trim_y(Fractures[i].Verts_trim.size()),
            Verts_trim_z(Fractures[i].Verts_trim.size());
        for (size_t j = 0; j < Fractures[i].Verts_trim.size(); ++j)
        {
            Verts_trim_x[j] = Fractures[i].Verts_trim[j][0];
            Verts_trim_y[j] = Fractures[i].Verts_trim[j][1];
            Verts_trim_z[j] = Fractures[i].Verts_trim[j][2];
        }

        vector<double> Center = {Fractures[i].Center[0],
                                 Fractures[i].Center[1],
                                 Fractures[i].Center[2]};

        vector<double> Normal_vector = {Fractures[i].Normal_vector[0],
                                        Fractures[i].Normal_vector[1],
                                        Fractures[i].Normal_vector[2]};

        vector<double> Plane_parameter = {Fractures[i].Plane_parameter[0],
                                          Fractures[i].Plane_parameter[1],
                                          Fractures[i].Plane_parameter[2],
                                          Fractures[i].Plane_parameter[3]};

        //-----------------------------------
        vector<string> string_field = {Frac_name + "Tag",
                                       Frac_name + "Clus",
                                       Frac_name + "Nvertices",
                                       Frac_name + "Nvertices_trim",
                                       Frac_name + "Radius",
                                       Frac_name + "Dip_direction",
                                       Frac_name + "Dip_angle",
                                       Frac_name + "Area",
                                       Frac_name + "Perimeter",
                                       Frac_name + "Area_trim",
                                       Frac_name + "Perimeter_trim",
                                       Frac_name + "Conductivity",
                                       Frac_name + "Verts_x",
                                       Frac_name + "Verts_y",
                                       Frac_name + "Verts_z",
                                       Frac_name + "Verts_trim_x",
                                       Frac_name + "Verts_trim_y",
                                       Frac_name + "Verts_trim_z",
                                       Frac_name + "Center",
                                       Frac_name + "Normal_vector",
                                       Frac_name + "Plane_parameter"};

        vector<vector<double>> data_set = { vector<double>{(double)Fractures[i].Tag},
                                            vector<double>{(double)Fractures[i].Clus},
                                            vector<double>{(double)Fractures[i].Nvertices},
                                            vector<double>{(double)Fractures[i].Nvertices_trim},
                                            vector<double>{(double)Fractures[i].Radius},
                                            vector<double>{(double)Fractures[i].Dip_direction},
                                            vector<double>{(double)Fractures[i].Dip_angle},
                                            vector<double>{(double)Fractures[i].Area},
                                            vector<double>{(double)Fractures[i].Perimeter},
                                            vector<double>{(double)Fractures[i].Area_trim},
                                            vector<double>{(double)Fractures[i].Perimeter_trim},
                                            vector<double>{(double)Fractures[i].Conductivity},
                                            Verts_x,
                                            Verts_y,
                                            Verts_z,
                                            Verts_trim_x,
                                            Verts_trim_y,
                                            Verts_trim_z,
                                            Center,
                                            Normal_vector,
                                            Plane_parameter };
        string_field_overall.insert(string_field_overall.end(), string_field.begin(), string_field.end());
        data_set_overall.insert(data_set_overall.end(), data_set.begin(), data_set.end());
    }
    h5file.Write_H5(filename, groupname, string_field_overall, data_set_overall);
    h5file.Append_dataset_to_group(filename, groupname, "Frac_Num", vector<double>{(double)Fractures.size()});
};
}; // namespace DFN