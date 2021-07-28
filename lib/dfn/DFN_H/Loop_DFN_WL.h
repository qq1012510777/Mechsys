#pragma once
#include "../FEM_H/FEM_DFN_A.h"
#include "../Mesh_H/Mesh_DFN_overall.h"
#include "Domain_WL.h"
#include <omp.h>
#include <stdio.h>
#include <unistd.h>

namespace DFN
{
class Loop_DFN
{
public:
    double times;                   ///< loop times, each time DFN modeling, the density will be increased compared to last time DFN modeling
    std::vector<Vector4d> array12;  ///< alpha (power law), min_radius, max_radius,
    double L;                       ///< model size
    size_t Nproc;                   ///< num of threads
    size_t nt;                      ///< should not larger than times! which model is the model shown to us
    size_t nk;                      ///< when nk DFN models are finished, output one time
    size_t nv_MC_TIMES;             ///< each density, the MC times
    double nx;                      ///< the increment of fracture number regard to each DFN modeling time
    double Percolation_parameter_c; ///< when percolation probability equals to 0.5, the percolation parameter is
    size_t NumofFsets;              ///< number of fracture sets
    size_t Nb_flow_sim_MC_times;

    std::vector<double> DenWeight; ///< the weight of number of each sets
    std::vector<Vector7d> array13; ///< the size of this array is the number of fracture sets;
    //and each element includes the seven inputs of each Fisher fracture set,
    //inputs are mean dip direction, mean dip angle, Fisher constant, min dip direction, max dip direction, min dip angle, max dip angle

    string Data_CommandFile = "Data_Command.m";
    string Data_MatFile = "Data_Mat.mat";
    /*
    std::vector<double> P32_total_1; ///< the next several 1D arrays store outputs
    std::vector<double> P32_connected_1;
    std::vector<double> P30_1;
    std::vector<double> Percolation_parameter_1;
    std::vector<double> Percolation_parameter_2;
    std::vector<double> Percolation_parameter_3;
    std::vector<double> Percolation_parameter_4;
    std::vector<double> Percolation_parameter_5;
    std::vector<double> Ratio_of_P32_1;
    std::vector<double> Percolation_probability_1;
    std::vector<double> n_I_1;
    std::vector<double> Correlation_length_1;
    std::vector<double> Max_gyration_radius_1;
    std::vector<double> P30_largest_cluster_1;
    std::vector<double> P32_largest_cluster_1;
    std::vector<double> P30_connected_1;
    std::vector<double> Ratio_of_P30_1;*/

public:
    void Loop_create_DFNs(gsl_rng *random_seed,
                          string str_ori,
                          string str_frac_size,
                          string percolation_direction,
                          const double min_ele_edge,
                          const double max_ele_edge,
                          const string conductivity_distri,
                          size_t modelno);

    void Data_output_stepBYstep(size_t times,
                                string FileKey,
                                double min_R,
                                double max_R,
                                double L,
                                double increment_fracture,
                                double P32_total_B,
                                double P32_connected_B,
                                double P30_B,
                                double Ratio_of_P32_B,
                                double Percolation_parameter_B_1,
                                double Percolation_parameter_B_2,
                                double Percolation_parameter_B_3,
                                double Percolation_parameter_B_4,
                                double Percolation_parameter_B_5,
                                double n_I,
                                double Percolation_probability_B,
                                double Correlation_length_B,
                                double Gyration_radius_B,
                                double P30_largest_cluster_B,
                                double P32_largest_cluster_B,
                                double P30_connected_B,
                                double Ratio_of_P30_B);

    void Sign_of_finding_pc(string FileKey);
    ///< if Pc is found, outputs a file

    void Matlab_Data_output_stepBYstep(const size_t np,
                                       string FileKey_mat,
                                       std::vector<double> P32_total_A,
                                       std::vector<double> P32_connected_A,
                                       std::vector<double> P30_A,
                                       std::vector<double> Percolation_parameter_A_1,
                                       std::vector<double> Percolation_parameter_A_2,
                                       std::vector<double> Percolation_parameter_A_3,
                                       std::vector<double> Percolation_parameter_A_4,
                                       std::vector<double> Percolation_parameter_A_5,
                                       std::vector<double> Ratio_of_P32_A,
                                       std::vector<double> Percolation_probability_A,
                                       std::vector<double> n_I_A,
                                       std::vector<double> Correlation_length_A,
                                       std::vector<double> Max_gyration_radius_A,
                                       std::vector<double> P30_largest_cluster_A,
                                       std::vector<double> P32_largest_cluster_A,
                                       std::vector<double> P30_connected_A,
                                       std::vector<double> Ratio_of_P30_A,
                                       std::vector<double> Permeability_A,
                                       const string str_ori,
                                       const string str_frac_size,
                                       const string conductivity_distri,
                                       const double domain_size);

    void Matlab_command(string FileKey_m, string FileKey_mat, size_t np, size_t ny, size_t model_no);

    bool If_conduct_flow_sim(size_t nf);
};

//****************************
inline void Loop_DFN::Loop_create_DFNs(gsl_rng *random_seed,
                                       string str_ori,
                                       string str_frac_size,
                                       string percolation_direction,
                                       const double min_ele_edge,
                                       const double max_ele_edge,
                                       const string conductivity_distri,
                                       size_t modelno)
{
    bool if_probability_1 = false;

    size_t nv = nv_MC_TIMES;
    //each density, the MC times

    std::vector<Vector2d> array11;
    //min x, max x;
    //min y, max y;
    //min z, max z;
    //(they are DOMAIN size, not MODEL size!!!)
    /*
        R_up = array12[0][2];
        for (size_t i = 0; i < NumofFsets; ++i)
        {
            if (R_up <= array12[i][2])
                R_up = array12[i][2];
        }
        */

    //if boundary effect is considered,
    //please uncomment the follow part
    double R_up = 0;
    /*
    if (str_frac_size == "powerlaw")
    {
        R_up = array12[0][2];
    }
    else if (str_frac_size == "lognormal")
    {
        R_up = array12[0][3];
    }
    else if (str_frac_size == "uniform")
    {
        R_up = array12[0][1];
    }
    else if (str_frac_size == "single")
    {
        R_up = array12[0][0];
    }
    else
    {
        throw Error_throw_pause("Error! Please define fracture size distribution!\n");
    }
    */
    array11.resize(3);
    array11[0][0] = -L * 0.5 - R_up;
    array11[0][1] = L * 0.5 + R_up;
    array11[1][0] = -L * 0.5 - R_up;
    array11[1][1] = L * 0.5 + R_up;
    array11[2][0] = -L * 0.5 - R_up;
    array11[2][1] = L * 0.5 + R_up;

    //   double array13[7] = {30, 20, 20, 10, 100.1, 0.5, 30.1};
    Vector6d model_size;
    model_size << -L * 0.5,
        L * 0.5,
        -L * 0.5,
        L * 0.5,
        -L * 0.5,
        L * 0.5; // MODEL size, not DOMAIN

    size_t np = 0;
    size_t njk = 0;

    while (np < times)
    {
        if (if_probability_1 == true)
        {
            nv = 0.5 * nv_MC_TIMES;
        }

        np++;
        size_t n = np * nx;

        //double nf = 0;

        std::vector<double> P32_total_A;
        std::vector<double> P32_connected_A;
        std::vector<double> P30_A;
        std::vector<double> Percolation_parameter_A_1;
        std::vector<double> Percolation_parameter_A_2;
        std::vector<double> Percolation_parameter_A_3;
        std::vector<double> Percolation_parameter_A_4;
        std::vector<double> Percolation_parameter_A_5;
        std::vector<double> Ratio_of_P32_A;
        std::vector<double> Percolation_probability_A;
        std::vector<double> n_I_A;
        std::vector<double> Correlation_length_A;
        std::vector<double> Max_gyration_radius_A;
        std::vector<double> P30_largest_cluster_A;
        std::vector<double> P32_largest_cluster_A;
        std::vector<double> P30_connected_A;
        std::vector<double> Ratio_of_P30_A;
        std::vector<double> Permeability_A;

        P32_total_A.resize(nv);
        P32_connected_A.resize(nv);
        P30_A.resize(nv);
        Percolation_parameter_A_1.resize(nv);
        Percolation_parameter_A_2.resize(nv);
        Percolation_parameter_A_3.resize(nv);
        Percolation_parameter_A_4.resize(nv);
        Percolation_parameter_A_5.resize(nv);
        Ratio_of_P32_A.resize(nv);
        Percolation_probability_A.resize(nv);
        n_I_A.resize(nv);
        Correlation_length_A.resize(nv);
        Max_gyration_radius_A.resize(nv);
        P30_largest_cluster_A.resize(nv);
        P32_largest_cluster_A.resize(nv);
        P30_connected_A.resize(nv);
        Ratio_of_P30_A.resize(nv);
        Permeability_A.resize(Nb_flow_sim_MC_times);

#pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i = 0; i < nv; i++)
        {

        Regenerate_dfn:;
            try
            {
                DFN::Domain dom;
                //dom.Fractures.clear();

                dom.Create_whole_model(n,
                                       DenWeight,
                                       random_seed,
                                       model_size,
                                       str_ori,
                                       str_frac_size,
                                       array11,
                                       array12,
                                       array13,
                                       conductivity_distri);
                ///uniform means oritation data
                //are generated uniformly, so,
                //actually, array13 is input but not used

                size_t z = dom.Identify_percolation_clusters(percolation_direction);
                if (str_ori == "uniform")
                    dom.Connectivity_uniform_orientation(percolation_direction);
                else if (str_ori == "fisher")
                    dom.Connectivity_fisher_orientation(percolation_direction);
                else
                {
                    throw Error_throw_pause("Error! Please define the orientation distribution!\n");
                }

                P32_total_A[i] = (dom.P32_total);
                P32_connected_A[i] = (dom.P32_connected);
                P30_A[i] = (dom.P30);
                Percolation_parameter_A_1[i] = (dom.Percolation_parameter_a);
                Percolation_parameter_A_2[i] = (dom.Percolation_parameter_b);
                Percolation_parameter_A_3[i] = (dom.Percolation_parameter_c);
                Percolation_parameter_A_4[i] = (dom.Percolation_parameter_d);
                Percolation_parameter_A_5[i] = (dom.Percolation_parameter_e);
                Ratio_of_P32_A[i] = (dom.Ratio_of_P32);
                n_I_A[i] = (dom.n_I);
                Correlation_length_A[i] = (dom.Xi);
                Max_gyration_radius_A[i] = (dom.max_R_s);
                P30_largest_cluster_A[i] = (dom.P30_largest_cluster);
                P32_largest_cluster_A[i] = (dom.P32_largest_cluster);
                P30_connected_A[i] = (dom.P30_connected);
                Ratio_of_P30_A[i] = (dom.Ratio_of_P30);

                /*
                if (i == nv - 1)
                {
                    cout << "Frac num: " << dom.Fractures.size() << endl;
                    cout << "V: " << (dom.Model_domain(0) - dom.Model_domain(1)) * (dom.Model_domain(3) - dom.Model_domain(2)) * (dom.Model_domain(5) - dom.Model_domain(4)) << endl;
                    cout << "P30: " << dom.Fractures.size() / (double)((dom.Model_domain(0) - dom.Model_domain(1)) * (dom.Model_domain(3) - dom.Model_domain(2)) * (dom.Model_domain(5) - dom.Model_domain(4))) << endl;
                    cout << "P30_c: " << dom.P30 << endl;
                }*/

                if (z == 1)
                    Percolation_probability_A[i] = 1;
                else
                    Percolation_probability_A[i] = 0;
                /*
#pragma omp critical
                {
                    if (Percolation_probability_A[i] == 1)
                        nf++;
                }
*/
                if (np == nt && i == nv - 1)
                {

                    dom.PlotMatlab_DFN("tdfn01_DFN.m");
                    dom.PlotMatlab_DFN_trim("tdfn01_DFN_trim.m");
                    dom.PlotMatlab_DFN_and_Intersection("tdfn01_DFN_and_Intersections.m");
                    dom.PlotMatlab_ORI_SCATTER("tdfn01_ORI_SCATTER.m");
                    //dom.PlotMatlab_Traces_on_Model_surfaces("tdfn01_Trace_on_surfaces.m");
                    dom.PlotMatlab_DFN_Highlight_Cluster("tdfn01_DFN_Highlight_Cluster.m");
                    dom.PLotMatlab_DFN_Cluster_along_a_direction("tdfn01_DFN_Z_clusters.m", "z");
                    //dom.PlotMatlab_Radius_and_Area_kstest("tdfn01_DFN_Fracture_Radius_and_Area.m");
                    //dom.PlotMatlab_Radius_and_Perimeter_kstest("tdfn01_DFN_Fracture_Radius_and_Perimeter.m");
                    //dom.DataFile_Radius_AreaAndPerimeter("tdfn01_DFN_Radius_AreaAndPerimeter.txt");
                    dom.Matlab_Out_Frac_matfile("Fractures.mat");
                }

                if (str_frac_size == "powerlaw")
                {
                    if (i == nv / 2 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; Alpha: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                    if (i == nv - 1 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; Alpha: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                }
                else if (str_frac_size == "lognormal")
                {
                    if (i == nv / 2 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; mean: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                    if (i == nv - 1 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; mean: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                }
                else if (str_frac_size == "uniform")
                {
                    if (i == nv / 2 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; lower: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                    if (i == nv - 1 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; lower: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                }
                else if (str_frac_size == "single")
                {
                    if (i == nv / 2 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; single_size: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                    if (i == nv - 1 && np % nk == 0)
                    {
                        using namespace std;
                        std::cout << "The Model NO." << np << " has been created! "
                                  << "Times: " << i << "; single_size: " << array12[0][0] << "; thread: " << omp_get_thread_num() << std::endl;
                    }
                }

                //-----------------------------------------------------------------------

                if (z == 1 && i < Nb_flow_sim_MC_times)
                {
                    /*
                    dom.Re_identify_intersection_considering_trimmed_frac();
                    size_t z2 = dom.Identify_percolation_clusters(percolation_direction);
                    //dom.PlotMatlab_DFN_and_Intersection("tdfn01_DFN_and_Intersections_II.m");
                    if (z2 != z)
                    {
                        throw Error_throw_ignore("Trimmed fractures did not form at least a percolation cluster!\n");
                    }

                    DFN::Mesh_DFN_overall mesh(dom, min_ele_edge, max_ele_edge);
                    DFN::FEM_DFN_A CC(mesh, dom);
                    if (np == nt && i == nv - 1)
                    {
                        mesh.Matlab_plot("mesh_DFN.mat", "mesh_DFN.m", dom);
                        CC.matlab_plot("FEM_DFN.mat", "FEM_DFN.m", dom, mesh, CC.F_overall);
                    }
                    Permeability_A[i] = CC.Permeability;
                    */
                    Permeability_A[i] = 0;
                }
                else if (z == 0 && i < Nb_flow_sim_MC_times)
                {
                    Permeability_A[i] = 0;
                }
            }
            catch (Error_throw_pause e)
            {
                cout << "\033[31mPause now! Because:\n";
                cout << e.msg << "\033[0m" << endl;
                exit(0);
            }
            catch (Error_throw_ignore e)
            {
                //cout << "\033[33mRegenerate a DFN! Because:\n";
                //cout << e.msg << "\033[0m" << endl;
                goto Regenerate_dfn;
            }
            catch (bad_alloc &e)
            {
                cout << "\033[33mRegenerate a DFN! Because:\n"
                     << e.what() << "\033[0m" << endl;
                goto Regenerate_dfn;
            }
        }

        this->Matlab_Data_output_stepBYstep(np,
                                            Data_MatFile,
                                            P32_total_A,
                                            P32_connected_A,
                                            P30_A,
                                            Percolation_parameter_A_1,
                                            Percolation_parameter_A_2,
                                            Percolation_parameter_A_3,
                                            Percolation_parameter_A_4,
                                            Percolation_parameter_A_5,
                                            Ratio_of_P32_A,
                                            Percolation_probability_A,
                                            n_I_A,
                                            Correlation_length_A,
                                            Max_gyration_radius_A,
                                            P30_largest_cluster_A,
                                            P32_largest_cluster_A,
                                            P30_connected_A,
                                            Ratio_of_P30_A,
                                            Permeability_A,
                                            str_ori,
                                            str_frac_size,
                                            conductivity_distri,
                                            this->L);

        double Percolation_parameter_B_1 = 0;
        for (size_t i = 0; i < Percolation_parameter_A_1.size(); ++i)
        {
            Percolation_parameter_B_1 += Percolation_parameter_A_1[i];
        }

        size_t nf = 0;
        if (njk == 0)
            for (size_t i = 0; i < Percolation_probability_A.size(); ++i)
                if (Percolation_probability_A[i] == 1)
                    nf++;

        if (njk == 0 && (double)nf / nv > 0.49999 && if_probability_1 == false)
        {
            njk++;
            Percolation_parameter_c = Percolation_parameter_B_1 / Percolation_parameter_A_1.size();
            std::cout << "\n**********Found Pc**********\n\n";
            cout << "reduce MC times!\n";
            if_probability_1 = true;
            Sign_of_finding_pc("Pc_Found.txt");
        }

        if (njk != 0 && Percolation_parameter_B_1 / Percolation_parameter_A_1.size() > 2 * Percolation_parameter_c)
        {
            std::cout << "\n**********Found two times Pc**********\n\n";
            break;
        }
    };

    this->Matlab_command(Data_CommandFile, Data_MatFile, np, np, modelno);
    std::cout << "Loop finished!\n";
};

inline void Loop_DFN::Data_output_stepBYstep(size_t times,
                                             string FileKey,
                                             double min_R,
                                             double max_R,
                                             double L,
                                             double increment_fracture,
                                             double P32_total_B,
                                             double P32_connected_B,
                                             double P30_B,
                                             double Ratio_of_P32_B,
                                             double Percolation_parameter_B_1,
                                             double Percolation_parameter_B_2,
                                             double Percolation_parameter_B_3,
                                             double Percolation_parameter_B_4,
                                             double Percolation_parameter_B_5,
                                             double n_I,
                                             double Percolation_probability_B,
                                             double Correlation_length_B,
                                             double Gyration_radius_B,
                                             double P30_largest_cluster_B,
                                             double P32_largest_cluster_B,
                                             double P30_connected_B,
                                             double Ratio_of_P30_B)
{
    if (times == 1)
    {
        //Writing data
        std::ofstream oss(FileKey, ios::out);
        oss << "Model_edge"
            << "\t" << L << "\n";
        oss << "Fracture_increment:"
            << "\t" << increment_fracture << "\n";
        oss << "min_R"
            << "\t" << min_R << "\n";
        oss << "max_R"
            << "\t" << max_R << "\n";
        oss << "P32_total"
            << "\t"
            << "P32_connected"
            << "\t"
            << "P30"
            << "\t"
            << "Ratio_of_P32"
            << "\t"
            << "Percolation_parameter_1"
            << "\t"
            << "Percolation_parameter_2"
            << "\t"
            << "Percolation_parameter_3"
            << "\t"
            << "Percolation_parameter_4"
            << "\t"
            << "Percolation_parameter_5"
            << "\t"
            << "n_I"
            << "\t"
            << "Percolation_probability"
            << "\t"
            << "Correlation_length_over_L"
            << "\t"
            << "Max_gyration_radius_over_L"
            << "\t"
            << "P30_largest_cluster"
            << "\t"
            << "P32_largest_cluster"
            << "\t"
            << "P30_connected"
            << "\t"
            << "Ratio_of_P30"
            << "\n";
        oss << P32_total_B << "\t" << P32_connected_B << "\t" << P30_B << "\t" << Ratio_of_P32_B << "\t" << Percolation_parameter_B_1 << "\t" << Percolation_parameter_B_2 << "\t" << Percolation_parameter_B_3 << "\t" << Percolation_parameter_B_4 << "\t" << Percolation_parameter_B_5 << "\t" << n_I << "\t" << Percolation_probability_B << "\t" << Correlation_length_B / L << "\t" << Gyration_radius_B / L << "\t" << P30_largest_cluster_B << "\t" << P32_largest_cluster_B << "\t" << P30_connected_B << "\t" << Ratio_of_P30_B << "\n";

        oss.close();
    }
    else
    {
        //Writing data
        std::ofstream oss(FileKey, ios::app);
        oss << P32_total_B << "\t" << P32_connected_B << "\t" << P30_B << "\t" << Ratio_of_P32_B << "\t" << Percolation_parameter_B_1 << "\t" << Percolation_parameter_B_2 << "\t" << Percolation_parameter_B_3 << "\t" << Percolation_parameter_B_4 << "\t" << Percolation_parameter_B_5 << "\t" << n_I << "\t" << Percolation_probability_B << "\t" << Correlation_length_B / L << "\t" << Gyration_radius_B / L << "\t" << P30_largest_cluster_B << "\t" << P32_largest_cluster_B << "\t" << P30_connected_B << "\t" << Ratio_of_P30_B << "\n";

        oss.close();
    }
};

inline void Loop_DFN::Matlab_Data_output_stepBYstep(const size_t np,
                                                    string FileKey_mat,
                                                    std::vector<double> P32_total_A,
                                                    std::vector<double> P32_connected_A,
                                                    std::vector<double> P30_A,
                                                    std::vector<double> Percolation_parameter_A_1,
                                                    std::vector<double> Percolation_parameter_A_2,
                                                    std::vector<double> Percolation_parameter_A_3,
                                                    std::vector<double> Percolation_parameter_A_4,
                                                    std::vector<double> Percolation_parameter_A_5,
                                                    std::vector<double> Ratio_of_P32_A,
                                                    std::vector<double> Percolation_probability_A,
                                                    std::vector<double> n_I_A,
                                                    std::vector<double> Correlation_length_A,
                                                    std::vector<double> Max_gyration_radius_A,
                                                    std::vector<double> P30_largest_cluster_A,
                                                    std::vector<double> P32_largest_cluster_A,
                                                    std::vector<double> P30_connected_A,
                                                    std::vector<double> Ratio_of_P30_A,
                                                    std::vector<double> Permeability_A,
                                                    const string str_ori,
                                                    const string str_frac_size,
                                                    const string conductivity_distri,
                                                    const double domain_size)

{

    if (np == 1)
    {
        const char *filename = FileKey_mat.c_str();
        MATFile *pMatFile;
        pMatFile = matOpen(filename, "w");

        if (!pMatFile)
        {
            cout << "Loop times: " << np << endl;
            cout << "cannot create mat file in class Loop_DFN_WL\n";
            throw Error_throw_pause("cannot create mat file in class Loop_DFN_WL\n");
        }

        //---------------
        double *pData19, *pData20, *pData21, *pData22;
        mxArray *pMxArray19, *pMxArray20, *pMxArray21, *pMxArray22;

        pData19 = (double *)mxCalloc(1, sizeof(double));
        pMxArray19 = mxCreateDoubleMatrix(1, 1, mxREAL);

        if (str_ori == "uniform")
        {
            pData19[0] = 11;
        }
        else if (str_ori == "fisher")
        {
            pData19 = (double *)mxCalloc(this->array13.size() * 7, sizeof(double));
            pMxArray19 = mxCreateDoubleMatrix(this->array13.size(), 7, mxREAL);

            for (size_t j = 0; j < this->array13.size() * 7; ++j)
            {
                size_t k, l;
                k = ceil(j / this->array13.size()); // column
                l = j % this->array13.size();       // row

                pData19[j] = this->array13[l](k);
            }
        }

        pData20 = (double *)mxCalloc(this->array12.size() * 4, sizeof(double));
        pMxArray20 = mxCreateDoubleMatrix(this->array12.size(), 4, mxREAL);
        for (size_t j = 0; j < this->array12.size() * 4; ++j)
        {
            size_t k, l;
            k = ceil(j / this->array12.size()); // column
            l = j % this->array12.size();       // row

            pData20[j] = this->array12[l](k);
        }

        pData21 = (double *)mxCalloc(1, sizeof(double));
        pMxArray21 = mxCreateDoubleMatrix(1, 1, mxREAL);
        pData21[0] = 15;

        pData22 = (double *)mxCalloc(1, sizeof(double));
        pMxArray22 = mxCreateDoubleMatrix(1, 1, mxREAL);
        pData22[0] = domain_size;

        if (!pMxArray19 || !pMxArray20 || !pMxArray21 || !pMxArray22)
        {
            cout << "Loop times: " << np << endl;
            cout << "cannot create pMxArray in class Loop_DFN_WL\n"
                 << endl;
            throw Error_throw_pause("cannot create pMxArray in class Loop_DFN_WL\n");
        }
        if (!pData19 || !pData20 || !pData21 || !pData22)
        {
            cout << "Loop times: " << np << endl;
            cout << "cannot create pData in class Loop_DFN_WL\n";
            cout << pData19[0] << "\n";
            cout << pData20[0] << "\n";
            cout << pData21[0] << "\n";
            cout << pData22[0] << "\n";
            throw Error_throw_pause("cannot create pData in class Loop_DFN_WL\n");
        }

        const char *char_ori = str_ori.c_str();
        const char *char_frac_size = str_frac_size.c_str();

        string conductivity_str = conductivity_distri + "_conductivity_dis";
        const char *char_conductivity_distri = conductivity_str.c_str();

        const char *char_domain_size = "Domain_size";

        mxSetData(pMxArray19, pData19);
        mxSetData(pMxArray20, pData20);
        mxSetData(pMxArray21, pData21);
        mxSetData(pMxArray22, pData22);

        matPutVariable(pMatFile, char_ori, pMxArray19);
        matPutVariable(pMatFile, char_frac_size, pMxArray20);
        matPutVariable(pMatFile, char_conductivity_distri, pMxArray21);
        matPutVariable(pMatFile, char_domain_size, pMxArray22);

        mxFree(pData19);
        mxFree(pData20);
        mxFree(pData21);
        mxFree(pData22);
        //---------------

        double *pData1,
            *pData2,
            *pData3,
            *pData4,
            *pData5,
            *pData6,
            *pData7,
            *pData8,
            *pData9,
            *pData10,
            *pData11,
            *pData12,
            *pData13,
            *pData14,
            *pData15,
            *pData16,
            *pData17,
            *pData18;

        pData1 = (double *)mxCalloc(P32_total_A.size(), sizeof(double));
        pData2 = (double *)mxCalloc(P32_connected_A.size(), sizeof(double));
        pData3 = (double *)mxCalloc(P30_A.size(), sizeof(double));
        pData4 = (double *)mxCalloc(Percolation_parameter_A_1.size(), sizeof(double));
        pData5 = (double *)mxCalloc(Percolation_parameter_A_2.size(), sizeof(double));
        pData6 = (double *)mxCalloc(Percolation_parameter_A_3.size(), sizeof(double));
        pData7 = (double *)mxCalloc(Percolation_parameter_A_4.size(), sizeof(double));
        pData8 = (double *)mxCalloc(Percolation_parameter_A_5.size(), sizeof(double));
        pData9 = (double *)mxCalloc(Ratio_of_P32_A.size(), sizeof(double));
        pData10 = (double *)mxCalloc(Percolation_probability_A.size(), sizeof(double));
        pData11 = (double *)mxCalloc(n_I_A.size(), sizeof(double));
        pData12 = (double *)mxCalloc(Correlation_length_A.size(), sizeof(double));
        pData13 = (double *)mxCalloc(Max_gyration_radius_A.size(), sizeof(double));
        pData14 = (double *)mxCalloc(P30_largest_cluster_A.size(), sizeof(double));
        pData15 = (double *)mxCalloc(P32_largest_cluster_A.size(), sizeof(double));
        pData16 = (double *)mxCalloc(P30_connected_A.size(), sizeof(double));
        pData17 = (double *)mxCalloc(Ratio_of_P30_A.size(), sizeof(double));
        pData18 = (double *)mxCalloc(Permeability_A.size(), sizeof(double));

        mxArray *pMxArray1;
        mxArray *pMxArray2;
        mxArray *pMxArray3;
        mxArray *pMxArray4;
        mxArray *pMxArray5;
        mxArray *pMxArray6;
        mxArray *pMxArray7;
        mxArray *pMxArray8;
        mxArray *pMxArray9;
        mxArray *pMxArray10;
        mxArray *pMxArray11;
        mxArray *pMxArray12;
        mxArray *pMxArray13;
        mxArray *pMxArray14;
        mxArray *pMxArray15;
        mxArray *pMxArray16;
        mxArray *pMxArray17;
        mxArray *pMxArray18;

        pMxArray1 = mxCreateDoubleMatrix(P32_total_A.size(), 1, mxREAL);
        pMxArray2 = mxCreateDoubleMatrix(P32_connected_A.size(), 1, mxREAL);
        pMxArray3 = mxCreateDoubleMatrix(P30_A.size(), 1, mxREAL);
        pMxArray4 = mxCreateDoubleMatrix(Percolation_parameter_A_1.size(), 1, mxREAL);
        pMxArray5 = mxCreateDoubleMatrix(Percolation_parameter_A_2.size(), 1, mxREAL);
        pMxArray6 = mxCreateDoubleMatrix(Percolation_parameter_A_3.size(), 1, mxREAL);
        pMxArray7 = mxCreateDoubleMatrix(Percolation_parameter_A_4.size(), 1, mxREAL);
        pMxArray8 = mxCreateDoubleMatrix(Percolation_parameter_A_5.size(), 1, mxREAL);
        pMxArray9 = mxCreateDoubleMatrix(Ratio_of_P32_A.size(), 1, mxREAL);
        pMxArray10 = mxCreateDoubleMatrix(Percolation_probability_A.size(), 1, mxREAL);
        pMxArray11 = mxCreateDoubleMatrix(n_I_A.size(), 1, mxREAL);
        pMxArray12 = mxCreateDoubleMatrix(Correlation_length_A.size(), 1, mxREAL);
        pMxArray13 = mxCreateDoubleMatrix(Max_gyration_radius_A.size(), 1, mxREAL);
        pMxArray14 = mxCreateDoubleMatrix(P30_largest_cluster_A.size(), 1, mxREAL);
        pMxArray15 = mxCreateDoubleMatrix(P32_largest_cluster_A.size(), 1, mxREAL);
        pMxArray16 = mxCreateDoubleMatrix(P30_connected_A.size(), 1, mxREAL);
        pMxArray17 = mxCreateDoubleMatrix(Ratio_of_P30_A.size(), 1, mxREAL);
        pMxArray18 = mxCreateDoubleMatrix(Permeability_A.size(), 1, mxREAL);

        if (!pMxArray1 ||
            !pMxArray2 ||
            !pMxArray3 ||
            !pMxArray4 ||
            !pMxArray5 ||
            !pMxArray6 ||
            !pMxArray7 ||
            !pMxArray8 ||
            !pMxArray9 ||
            !pMxArray10 ||
            !pMxArray11 ||
            !pMxArray12 ||
            !pMxArray13 ||
            !pMxArray14 ||
            !pMxArray15 ||
            !pMxArray16 ||
            !pMxArray17 ||
            !pMxArray18)
        {
            cout << "Loop times: " << np << endl;
            throw Error_throw_pause("cannot create pMxArray in class Loop_DFN_WL\n");
        }

        if (!pData1 ||
            !pData2 ||
            !pData3 ||
            !pData4 ||
            !pData5 ||
            !pData6 ||
            !pData7 ||
            !pData8 ||
            !pData9 ||
            !pData10 ||
            !pData11 ||
            !pData12 ||
            !pData13 ||
            !pData14 ||
            !pData15 ||
            !pData16 ||
            !pData17 ||
            !pData18)
        {
            cout << "Loop times: " << np << endl;
            throw Error_throw_pause("cannot create pData in class Loop_DFN_WL\n");
        }

        for (size_t i = 0; i < P32_total_A.size(); ++i)
        {
            pData1[i] = P32_total_A[i];
            pData2[i] = P32_connected_A[i];
            pData3[i] = P30_A[i];
            pData4[i] = Percolation_parameter_A_1[i];
            pData5[i] = Percolation_parameter_A_2[i];
            pData6[i] = Percolation_parameter_A_3[i];
            pData7[i] = Percolation_parameter_A_4[i];
            pData8[i] = Percolation_parameter_A_5[i];
            pData9[i] = Ratio_of_P32_A[i];
            pData10[i] = Percolation_probability_A[i];
            pData11[i] = n_I_A[i];
            pData12[i] = Correlation_length_A[i];
            pData13[i] = Max_gyration_radius_A[i];
            pData14[i] = P30_largest_cluster_A[i];
            pData15[i] = P32_largest_cluster_A[i];
            pData16[i] = P30_connected_A[i];
            pData17[i] = Ratio_of_P30_A[i];
            if (i < Permeability_A.size())
                pData18[i] = Permeability_A[i];
        }

        mxSetData(pMxArray1, pData1);
        mxSetData(pMxArray2, pData2);
        mxSetData(pMxArray3, pData3);
        mxSetData(pMxArray4, pData4);
        mxSetData(pMxArray5, pData5);
        mxSetData(pMxArray6, pData6);
        mxSetData(pMxArray7, pData7);
        mxSetData(pMxArray8, pData8);
        mxSetData(pMxArray9, pData9);
        mxSetData(pMxArray10, pData10);
        mxSetData(pMxArray11, pData11);
        mxSetData(pMxArray12, pData12);
        mxSetData(pMxArray13, pData13);
        mxSetData(pMxArray14, pData14);
        mxSetData(pMxArray15, pData15);
        mxSetData(pMxArray16, pData16);
        mxSetData(pMxArray17, pData17);
        mxSetData(pMxArray18, pData18);

        string ft = to_string(np);

        string string_P32_total = "P32_total_" + ft;
        string string_P32_connected = "P32_connected_" + ft;
        string string_P30 = "P30_" + ft;
        string string_Percolation_parameter_A_1 = "Percolation_parameter_1A_" + ft;
        string string_Percolation_parameter_A_2 = "Percolation_parameter_2A_" + ft;
        string string_Percolation_parameter_A_3 = "Percolation_parameter_3A_" + ft;
        string string_Percolation_parameter_A_4 = "Percolation_parameter_4A_" + ft;
        string string_Percolation_parameter_A_5 = "Percolation_parameter_5A_" + ft;
        string string_Ratio_of_P32 = "Ratio_of_P32_" + ft;
        string string_Percolation_probability = "Percolation_probability_" + ft;
        string string_n_I = "n_I_" + ft;
        string string_Correlation_length = "Correlation_length_" + ft;
        string string_Max_gyration_radius = "Max_gyration_radius_" + ft;
        string string_P30_largest_cluster = "P30_largest_cluster_" + ft;
        string string_P32_largest_cluster = "P32_largest_cluster_" + ft;
        string string_P30_connected = "P30_connected_" + ft;
        string string_Ratio_of_P30 = "Ratio_of_P30_" + ft;
        string string_Permeability = "Permeability_" + ft;

        const char *char_P32_total = string_P32_total.c_str();
        const char *char_P32_connected = string_P32_connected.c_str();
        const char *char_P30 = string_P30.c_str();
        const char *char_Percolation_parameter_A_1 = string_Percolation_parameter_A_1.c_str();
        const char *char_Percolation_parameter_A_2 = string_Percolation_parameter_A_2.c_str();
        const char *char_Percolation_parameter_A_3 = string_Percolation_parameter_A_3.c_str();
        const char *char_Percolation_parameter_A_4 = string_Percolation_parameter_A_4.c_str();
        const char *char_Percolation_parameter_A_5 = string_Percolation_parameter_A_5.c_str();
        const char *char_Ratio_of_P32 = string_Ratio_of_P32.c_str();
        const char *char_Percolation_probability = string_Percolation_probability.c_str();
        const char *char_n_I = string_n_I.c_str();
        const char *char_Correlation_length = string_Correlation_length.c_str();
        const char *char_Max_gyration_radius = string_Max_gyration_radius.c_str();
        const char *char_P30_largest_cluster = string_P30_largest_cluster.c_str();
        const char *char_P32_largest_cluster = string_P32_largest_cluster.c_str();
        const char *char_P30_connected = string_P30_connected.c_str();
        const char *char_Ratio_of_P30 = string_Ratio_of_P30.c_str();
        const char *char_Permeability = string_Permeability.c_str();

        matPutVariable(pMatFile, char_P32_total, pMxArray1);
        matPutVariable(pMatFile, char_P32_connected, pMxArray2);
        matPutVariable(pMatFile, char_P30, pMxArray3);
        matPutVariable(pMatFile, char_Percolation_parameter_A_1, pMxArray4);
        matPutVariable(pMatFile, char_Percolation_parameter_A_2, pMxArray5);
        matPutVariable(pMatFile, char_Percolation_parameter_A_3, pMxArray6);
        matPutVariable(pMatFile, char_Percolation_parameter_A_4, pMxArray7);
        matPutVariable(pMatFile, char_Percolation_parameter_A_5, pMxArray8);
        matPutVariable(pMatFile, char_Ratio_of_P32, pMxArray9);
        matPutVariable(pMatFile, char_Percolation_probability, pMxArray10);
        matPutVariable(pMatFile, char_n_I, pMxArray11);
        matPutVariable(pMatFile, char_Correlation_length, pMxArray12);
        matPutVariable(pMatFile, char_Max_gyration_radius, pMxArray13);
        matPutVariable(pMatFile, char_P30_largest_cluster, pMxArray14);
        matPutVariable(pMatFile, char_P32_largest_cluster, pMxArray15);
        matPutVariable(pMatFile, char_P30_connected, pMxArray16);
        matPutVariable(pMatFile, char_Ratio_of_P30, pMxArray17);
        matPutVariable(pMatFile, char_Permeability, pMxArray18);

        mxFree(pData1);
        mxFree(pData2);
        mxFree(pData3);
        mxFree(pData4);
        mxFree(pData5);
        mxFree(pData6);
        mxFree(pData7);
        mxFree(pData8);
        mxFree(pData9);
        mxFree(pData10);
        mxFree(pData11);
        mxFree(pData12);
        mxFree(pData13);
        mxFree(pData14);
        mxFree(pData15);
        mxFree(pData16);
        mxFree(pData17);

        matClose(pMatFile);
    }
    else
    {
        const char *filename = FileKey_mat.c_str();
        MATFile *pMatFile;
        pMatFile = matOpen(filename, "u");

        if (!pMatFile)
        {
            cout << "Loop times: " << np << endl;
            throw Error_throw_pause("cannot create mat file in class Loop_DFN_WL\n");
        }

        double *pData1,
            *pData2,
            *pData3,
            *pData4,
            *pData5,
            *pData6,
            *pData7,
            *pData8,
            *pData9,
            *pData10,
            *pData11,
            *pData12,
            *pData13,
            *pData14,
            *pData15,
            *pData16,
            *pData17,
            *pData18;

        pData1 = (double *)mxCalloc(P32_total_A.size(), sizeof(double));
        pData2 = (double *)mxCalloc(P32_connected_A.size(), sizeof(double));
        pData3 = (double *)mxCalloc(P30_A.size(), sizeof(double));
        pData4 = (double *)mxCalloc(Percolation_parameter_A_1.size(), sizeof(double));
        pData5 = (double *)mxCalloc(Percolation_parameter_A_2.size(), sizeof(double));
        pData6 = (double *)mxCalloc(Percolation_parameter_A_3.size(), sizeof(double));
        pData7 = (double *)mxCalloc(Percolation_parameter_A_4.size(), sizeof(double));
        pData8 = (double *)mxCalloc(Percolation_parameter_A_5.size(), sizeof(double));
        pData9 = (double *)mxCalloc(Ratio_of_P32_A.size(), sizeof(double));
        pData10 = (double *)mxCalloc(Percolation_probability_A.size(), sizeof(double));
        pData11 = (double *)mxCalloc(n_I_A.size(), sizeof(double));
        pData12 = (double *)mxCalloc(Correlation_length_A.size(), sizeof(double));
        pData13 = (double *)mxCalloc(Max_gyration_radius_A.size(), sizeof(double));
        pData14 = (double *)mxCalloc(P30_largest_cluster_A.size(), sizeof(double));
        pData15 = (double *)mxCalloc(P32_largest_cluster_A.size(), sizeof(double));
        pData16 = (double *)mxCalloc(P30_connected_A.size(), sizeof(double));
        pData17 = (double *)mxCalloc(Ratio_of_P30_A.size(), sizeof(double));
        pData18 = (double *)mxCalloc(Permeability_A.size(), sizeof(double));

        mxArray *pMxArray1;
        mxArray *pMxArray2;
        mxArray *pMxArray3;
        mxArray *pMxArray4;
        mxArray *pMxArray5;
        mxArray *pMxArray6;
        mxArray *pMxArray7;
        mxArray *pMxArray8;
        mxArray *pMxArray9;
        mxArray *pMxArray10;
        mxArray *pMxArray11;
        mxArray *pMxArray12;
        mxArray *pMxArray13;
        mxArray *pMxArray14;
        mxArray *pMxArray15;
        mxArray *pMxArray16;
        mxArray *pMxArray17;
        mxArray *pMxArray18;

        pMxArray1 = mxCreateDoubleMatrix(P32_total_A.size(), 1, mxREAL);
        pMxArray2 = mxCreateDoubleMatrix(P32_connected_A.size(), 1, mxREAL);
        pMxArray3 = mxCreateDoubleMatrix(P30_A.size(), 1, mxREAL);
        pMxArray4 = mxCreateDoubleMatrix(Percolation_parameter_A_1.size(), 1, mxREAL);
        pMxArray5 = mxCreateDoubleMatrix(Percolation_parameter_A_2.size(), 1, mxREAL);
        pMxArray6 = mxCreateDoubleMatrix(Percolation_parameter_A_3.size(), 1, mxREAL);
        pMxArray7 = mxCreateDoubleMatrix(Percolation_parameter_A_4.size(), 1, mxREAL);
        pMxArray8 = mxCreateDoubleMatrix(Percolation_parameter_A_5.size(), 1, mxREAL);
        pMxArray9 = mxCreateDoubleMatrix(Ratio_of_P32_A.size(), 1, mxREAL);
        pMxArray10 = mxCreateDoubleMatrix(Percolation_probability_A.size(), 1, mxREAL);
        pMxArray11 = mxCreateDoubleMatrix(n_I_A.size(), 1, mxREAL);
        pMxArray12 = mxCreateDoubleMatrix(Correlation_length_A.size(), 1, mxREAL);
        pMxArray13 = mxCreateDoubleMatrix(Max_gyration_radius_A.size(), 1, mxREAL);
        pMxArray14 = mxCreateDoubleMatrix(P30_largest_cluster_A.size(), 1, mxREAL);
        pMxArray15 = mxCreateDoubleMatrix(P32_largest_cluster_A.size(), 1, mxREAL);
        pMxArray16 = mxCreateDoubleMatrix(P30_connected_A.size(), 1, mxREAL);
        pMxArray17 = mxCreateDoubleMatrix(Ratio_of_P30_A.size(), 1, mxREAL);
        pMxArray18 = mxCreateDoubleMatrix(Permeability_A.size(), 1, mxREAL);

        if (!pMxArray1 ||
            !pMxArray2 ||
            !pMxArray3 ||
            !pMxArray4 ||
            !pMxArray5 ||
            !pMxArray6 ||
            !pMxArray7 ||
            !pMxArray8 ||
            !pMxArray9 ||
            !pMxArray10 ||
            !pMxArray11 ||
            !pMxArray12 ||
            !pMxArray13 ||
            !pMxArray14 ||
            !pMxArray15 ||
            !pMxArray16 ||
            !pMxArray17 ||
            !pMxArray18)
        {
            cout << "Loop times: " << np << endl;
            throw Error_throw_pause("cannot create pMxArray in class Loop_DFN_WL\n");
        }

        if (!pData1 ||
            !pData2 ||
            !pData3 ||
            !pData4 ||
            !pData5 ||
            !pData6 ||
            !pData7 ||
            !pData8 ||
            !pData9 ||
            !pData10 ||
            !pData11 ||
            !pData12 ||
            !pData13 ||
            !pData14 ||
            !pData15 ||
            !pData16 ||
            !pData17 ||
            !pData18)
        {
            cout << "Loop times: " << np << endl;
            throw Error_throw_pause("cannot create pData in class Loop_DFN_WL\n");
        }

        for (size_t i = 0; i < P32_total_A.size(); ++i)
        {
            pData1[i] = P32_total_A[i];
            pData2[i] = P32_connected_A[i];
            pData3[i] = P30_A[i];
            pData4[i] = Percolation_parameter_A_1[i];
            pData5[i] = Percolation_parameter_A_2[i];
            pData6[i] = Percolation_parameter_A_3[i];
            pData7[i] = Percolation_parameter_A_4[i];
            pData8[i] = Percolation_parameter_A_5[i];
            pData9[i] = Ratio_of_P32_A[i];
            pData10[i] = Percolation_probability_A[i];
            pData11[i] = n_I_A[i];
            pData12[i] = Correlation_length_A[i];
            pData13[i] = Max_gyration_radius_A[i];
            pData14[i] = P30_largest_cluster_A[i];
            pData15[i] = P32_largest_cluster_A[i];
            pData16[i] = P30_connected_A[i];
            pData17[i] = Ratio_of_P30_A[i];
            if (i < Permeability_A.size())
                pData18[i] = Permeability_A[i];
        }
        mxSetData(pMxArray1, pData1);
        mxSetData(pMxArray2, pData2);
        mxSetData(pMxArray3, pData3);
        mxSetData(pMxArray4, pData4);
        mxSetData(pMxArray5, pData5);
        mxSetData(pMxArray6, pData6);
        mxSetData(pMxArray7, pData7);
        mxSetData(pMxArray8, pData8);
        mxSetData(pMxArray9, pData9);
        mxSetData(pMxArray10, pData10);
        mxSetData(pMxArray11, pData11);
        mxSetData(pMxArray12, pData12);
        mxSetData(pMxArray13, pData13);
        mxSetData(pMxArray14, pData14);
        mxSetData(pMxArray15, pData15);
        mxSetData(pMxArray16, pData16);
        mxSetData(pMxArray17, pData17);
        mxSetData(pMxArray18, pData18);

        string ft = to_string(np);

        string string_P32_total = "P32_total_" + ft;
        string string_P32_connected = "P32_connected_" + ft;
        string string_P30 = "P30_" + ft;
        string string_Percolation_parameter_A_1 = "Percolation_parameter_1A_" + ft;
        string string_Percolation_parameter_A_2 = "Percolation_parameter_2A_" + ft;
        string string_Percolation_parameter_A_3 = "Percolation_parameter_3A_" + ft;
        string string_Percolation_parameter_A_4 = "Percolation_parameter_4A_" + ft;
        string string_Percolation_parameter_A_5 = "Percolation_parameter_5A_" + ft;
        string string_Ratio_of_P32 = "Ratio_of_P32_" + ft;
        string string_Percolation_probability = "Percolation_probability_" + ft;
        string string_n_I = "n_I_" + ft;
        string string_Correlation_length = "Correlation_length_" + ft;
        string string_Max_gyration_radius = "Max_gyration_radius_" + ft;
        string string_P30_largest_cluster = "P30_largest_cluster_" + ft;
        string string_P32_largest_cluster = "P32_largest_cluster_" + ft;
        string string_P30_connected = "P30_connected_" + ft;
        string string_Ratio_of_P30 = "Ratio_of_P30_" + ft;
        string string_Permeability = "Permeability_" + ft;

        const char *char_P32_total = string_P32_total.c_str();
        const char *char_P32_connected = string_P32_connected.c_str();
        const char *char_P30 = string_P30.c_str();
        const char *char_Percolation_parameter_A_1 = string_Percolation_parameter_A_1.c_str();
        const char *char_Percolation_parameter_A_2 = string_Percolation_parameter_A_2.c_str();
        const char *char_Percolation_parameter_A_3 = string_Percolation_parameter_A_3.c_str();
        const char *char_Percolation_parameter_A_4 = string_Percolation_parameter_A_4.c_str();
        const char *char_Percolation_parameter_A_5 = string_Percolation_parameter_A_5.c_str();
        const char *char_Ratio_of_P32 = string_Ratio_of_P32.c_str();
        const char *char_Percolation_probability = string_Percolation_probability.c_str();
        const char *char_n_I = string_n_I.c_str();
        const char *char_Correlation_length = string_Correlation_length.c_str();
        const char *char_Max_gyration_radius = string_Max_gyration_radius.c_str();
        const char *char_P30_largest_cluster = string_P30_largest_cluster.c_str();
        const char *char_P32_largest_cluster = string_P32_largest_cluster.c_str();
        const char *char_P30_connected = string_P30_connected.c_str();
        const char *char_Ratio_of_P30 = string_Ratio_of_P30.c_str();
        const char *char_Permeability = string_Permeability.c_str();

        matPutVariable(pMatFile, char_P32_total, pMxArray1);
        matPutVariable(pMatFile, char_P32_connected, pMxArray2);
        matPutVariable(pMatFile, char_P30, pMxArray3);
        matPutVariable(pMatFile, char_Percolation_parameter_A_1, pMxArray4);
        matPutVariable(pMatFile, char_Percolation_parameter_A_2, pMxArray5);
        matPutVariable(pMatFile, char_Percolation_parameter_A_3, pMxArray6);
        matPutVariable(pMatFile, char_Percolation_parameter_A_4, pMxArray7);
        matPutVariable(pMatFile, char_Percolation_parameter_A_5, pMxArray8);
        matPutVariable(pMatFile, char_Ratio_of_P32, pMxArray9);
        matPutVariable(pMatFile, char_Percolation_probability, pMxArray10);
        matPutVariable(pMatFile, char_n_I, pMxArray11);
        matPutVariable(pMatFile, char_Correlation_length, pMxArray12);
        matPutVariable(pMatFile, char_Max_gyration_radius, pMxArray13);
        matPutVariable(pMatFile, char_P30_largest_cluster, pMxArray14);
        matPutVariable(pMatFile, char_P32_largest_cluster, pMxArray15);
        matPutVariable(pMatFile, char_P30_connected, pMxArray16);
        matPutVariable(pMatFile, char_Ratio_of_P30, pMxArray17);
        matPutVariable(pMatFile, char_Permeability, pMxArray18);

        mxFree(pData1);
        mxFree(pData2);
        mxFree(pData3);
        mxFree(pData4);
        mxFree(pData5);
        mxFree(pData6);
        mxFree(pData7);
        mxFree(pData8);
        mxFree(pData9);
        mxFree(pData10);
        mxFree(pData11);
        mxFree(pData12);
        mxFree(pData13);
        mxFree(pData14);
        mxFree(pData15);
        mxFree(pData16);
        mxFree(pData17);

        matClose(pMatFile);
    }
};

inline void Loop_DFN::Matlab_command(string FileKey_m, string FileKey_mat, size_t np, size_t ny, size_t model_no)
{
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "s_" << model_no << " = load('" << FileKey_mat << "');\n";

    string L = "L_" + to_string(model_no) + " = s_" + to_string(model_no) + ".Domain_size;\n";
    oss << L;
    string string_P32_total = "P32_total_" + to_string(model_no);
    string string_P30 = "P30_" + to_string(model_no);
    string string_P32_connected = "P32_connected_" + to_string(model_no);
    string string_Ratio_of_P32 = "Ratio_of_P32_" + to_string(model_no);
    string string_Percolation_parameter_A_1 = "Percolation_parameter_1_" + to_string(model_no);
    string string_Percolation_parameter_A_2 = "Percolation_parameter_2_" + to_string(model_no);
    string string_Percolation_parameter_A_3 = "Percolation_parameter_3_" + to_string(model_no);
    string string_Percolation_parameter_A_4 = "Percolation_parameter_4_" + to_string(model_no);
    string string_Percolation_parameter_A_5 = "Percolation_parameter_5_" + to_string(model_no);
    string string_n_I = "n_I_" + to_string(model_no);
    string string_Percolation_probability = "Percolation_probability_" + to_string(model_no);
    string string_Correlation_length = "Correlation_length_" + to_string(model_no);
    string string_Max_gyration_radius = "Max_gyration_radius_" + to_string(model_no);
    string string_P30_largest_cluster = "P30_largest_cluster_" + to_string(model_no);
    string string_P32_largest_cluster = "P32_largest_cluster_" + to_string(model_no);
    string string_P30_connected = "P30_connected_" + to_string(model_no);
    string string_Ratio_of_P30 = "Ratio_of_P30_" + to_string(model_no);
    string string_Permeability = "Permeability_" + to_string(model_no);

    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_P32_total << " = [mean(s_" << model_no << ".P32_total_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".P32_total_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".P32_total_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_P30 << " = [mean(s_" << model_no << ".P30_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".P30_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".P30_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_P32_connected << " = [mean(s_" << model_no << ".P32_connected_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".P32_connected_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".P32_connected_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Ratio_of_P32 << " = [mean(s_" << model_no << ".Ratio_of_P32_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Ratio_of_P32_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Ratio_of_P32_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Percolation_parameter_A_1 << " = [mean(s_" << model_no << ".Percolation_parameter_1A_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Percolation_parameter_1A_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Percolation_parameter_1A_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Percolation_parameter_A_2 << " = [mean(s_" << model_no << ".Percolation_parameter_2A_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Percolation_parameter_2A_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Percolation_parameter_2A_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Percolation_parameter_A_3 << " = [mean(s_" << model_no << ".Percolation_parameter_3A_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Percolation_parameter_3A_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Percolation_parameter_3A_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Percolation_parameter_A_4 << " = [mean(s_" << model_no << ".Percolation_parameter_4A_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Percolation_parameter_4A_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Percolation_parameter_4A_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Percolation_parameter_A_5 << " = [mean(s_" << model_no << ".Percolation_parameter_5A_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Percolation_parameter_5A_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Percolation_parameter_5A_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_n_I << " = [mean(s_" << model_no << ".n_I_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".n_I_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".n_I_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Percolation_probability << " = [mean(s_" << model_no << ".Percolation_probability_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Percolation_probability_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Percolation_probability_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Correlation_length << " = [mean(s_" << model_no << ".Correlation_length_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Correlation_length_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Correlation_length_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Max_gyration_radius << " = [mean(s_" << model_no << ".Max_gyration_radius_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Max_gyration_radius_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Max_gyration_radius_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_P30_largest_cluster << " = [mean(s_" << model_no << ".P30_largest_cluster_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".P30_largest_cluster_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".P30_largest_cluster_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_P32_largest_cluster << " = [mean(s_" << model_no << ".P32_largest_cluster_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".P32_largest_cluster_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".P32_largest_cluster_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_P30_connected << " = [mean(s_" << model_no << ".P30_connected_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".P30_connected_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".P30_connected_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Ratio_of_P30 << " = [mean(s_" << model_no << ".Ratio_of_P30_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Ratio_of_P30_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Ratio_of_P30_" << ft << "); ...\n";
    }
    for (size_t ft = 1; ft <= np; ++ft)
    {
        if (ft == 1)
            oss << string_Permeability << " = [mean(s_" << model_no << ".Permeability_" << ft << "); ...\n";
        else if (ft == np)
            oss << "mean(s_" << model_no << ".Permeability_" << ft << ")];\n";
        else
            oss << "mean(s_" << model_no << ".Permeability_" << ft << "); ...\n";
    }

    oss << "clear s_" << model_no << ";\n";
    oss << "filepath = which('" << FileKey_m << "');\n";
    oss << "filepath = erase(filepath, '" << FileKey_m << "');\n";
    oss << "filepath = [filepath, 'connectivity_modelno_" << model_no << ".mat'];\n";
    oss << "save(filepath);\n";

    oss.close();
};

inline void Loop_DFN::Sign_of_finding_pc(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "Pc has been found: ";
    oss << Percolation_parameter_c;

    oss.close();
};

}; // namespace DFN
