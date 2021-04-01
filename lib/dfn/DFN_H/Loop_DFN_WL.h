#pragma once
#include "Domain_WL.h"
#include "FEM_DFN_WL.h"
#include "Mesh_DFN_WL.h"
#include <omp.h>

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

    std::vector<double> DenWeight;   ///< the weight of number of each sets
    std::vector<Vector7d> array13;   ///< the size of this array is the number of fracture sets; and each element includes the seven inputs of each Fisher fracture set, inputs are mean dip direction, mean dip angle, Fisher constant, min dip direction, max dip direction, min dip angle, max dip angle
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
    std::vector<double> Ratio_of_P30_1;

public:
    void Loop_create_DFNs(gsl_rng *random_seed,
                          string str_ori,
                          string str_frac_size,
                          string percolation_direction,
                          double avg_ele_len,
                          double ratio_H);

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
};

//****************************
inline void Loop_DFN::Loop_create_DFNs(gsl_rng *random_seed,
                                       string str_ori,
                                       string str_frac_size,
                                       string percolation_direction,
                                       double avg_ele_len,
                                       double ratio_H)
{
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
        std::cout << "Error! Please define fracture size distribution!\n";
        exit(0);
    }

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
        np++;
        size_t n = np * nx;

        double nf = 0;

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

#pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i = 0; i < nv; i++)
        {
            DFN::Domain dom;
            //std::cout<<"debug1\n";

            dom.Create_whole_model(n,
                                   DenWeight,
                                   random_seed,
                                   model_size,
                                   str_ori,
                                   str_frac_size,
                                   array11,
                                   array12,
                                   array13);
            ///uniform means oritation data
            //are generated uniformly, so,
            //actually, array13 is input but not used
            //std::cout<<"debug2\n";

            size_t z = dom.Identify_percolation_clusters(percolation_direction);
            if (str_ori == "uniform")
                dom.Connectivity_uniform_orientation(percolation_direction);
            else if (str_ori == "fisher")
                dom.Connectivity_fisher_orientation(percolation_direction);
            else
            {
                std::cout << "Error! Please define the orientation distribution!\n";
                exit(0);
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

            if (z == 1)
                Percolation_probability_A[i] = 1;
            else
                Percolation_probability_A[i] = 0;

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

            /*
            if (z == 1)
            {
                DFN::DFN_mesh AA(dom, avg_ele_len, ratio_H, percolation_direction);
                //std::cout << "mesh finished\n";
                DFN::FEM_DFN CC(AA, dom);
            }
            */
        }

        double P32_total_B = 0;
        for (size_t i = 0; i < P32_total_A.size(); ++i)
        {
            P32_total_B = P32_total_B + P32_total_A[i];
        }
        P32_total_1.push_back(P32_total_B / P32_total_A.size());

        double P32_connected_B = 0;
        for (size_t i = 0; i < P32_connected_A.size(); ++i)
        {
            P32_connected_B = P32_connected_B + P32_connected_A[i];
        }
        P32_connected_1.push_back(P32_connected_B / P32_connected_A.size());

        double P30_B = 0;
        for (size_t i = 0; i < P30_A.size(); ++i)
        {
            P30_B = P30_B + P30_A[i];
        }
        P30_1.push_back(P30_B / P30_A.size());

        double Percolation_parameter_B_1 = 0;
        for (size_t i = 0; i < Percolation_parameter_A_1.size(); ++i)
        {
            Percolation_parameter_B_1 = Percolation_parameter_B_1 + Percolation_parameter_A_1[i];
        }
        Percolation_parameter_1.push_back(Percolation_parameter_B_1 / Percolation_parameter_A_1.size());

        double Percolation_parameter_B_2 = 0;
        for (size_t i = 0; i < Percolation_parameter_A_2.size(); ++i)
        {
            Percolation_parameter_B_2 = Percolation_parameter_B_2 + Percolation_parameter_A_2[i];
        }
        Percolation_parameter_2.push_back(Percolation_parameter_B_2 / Percolation_parameter_A_2.size());

        double Percolation_parameter_B_3 = 0;
        for (size_t i = 0; i < Percolation_parameter_A_3.size(); ++i)
        {
            Percolation_parameter_B_3 = Percolation_parameter_B_3 + Percolation_parameter_A_3[i];
        }
        Percolation_parameter_3.push_back(Percolation_parameter_B_3 / Percolation_parameter_A_3.size());

        double Percolation_parameter_B_4 = 0;
        for (size_t i = 0; i < Percolation_parameter_A_4.size(); ++i)
        {
            Percolation_parameter_B_4 = Percolation_parameter_B_4 + Percolation_parameter_A_4[i];
        }
        Percolation_parameter_4.push_back(Percolation_parameter_B_4 / Percolation_parameter_A_4.size());

        double Percolation_parameter_B_5 = 0;
        for (size_t i = 0; i < Percolation_parameter_A_5.size(); ++i)
        {
            Percolation_parameter_B_5 = Percolation_parameter_B_5 + Percolation_parameter_A_5[i];
        }
        Percolation_parameter_5.push_back(Percolation_parameter_B_5 / Percolation_parameter_A_5.size());

        double Ratio_of_P32_B = 0;
        for (size_t i = 0; i < Ratio_of_P32_A.size(); ++i)
        {
            Ratio_of_P32_B = Ratio_of_P32_B + Ratio_of_P32_A[i];
        }
        Ratio_of_P32_1.push_back(Ratio_of_P32_B / Ratio_of_P32_A.size());

        double Correlation_length_B = 0;
        for (size_t i = 0; i < Correlation_length_A.size(); ++i)
        {
            Correlation_length_B = Correlation_length_B + Correlation_length_A[i];
        }
        Correlation_length_1.push_back(Correlation_length_B / Correlation_length_A.size());

        double Max_gyration_radius_B = 0;
        for (size_t i = 0; i < Max_gyration_radius_A.size(); ++i)
        {
            Max_gyration_radius_B = Max_gyration_radius_B + Max_gyration_radius_A[i];
        }
        Max_gyration_radius_1.push_back(Max_gyration_radius_B / Max_gyration_radius_A.size());

        double P30_largest_cluster_B = 0;
        for (size_t i = 0; i < P30_largest_cluster_A.size(); ++i)
        {
            P30_largest_cluster_B = P30_largest_cluster_B + P30_largest_cluster_A[i];
        }
        P30_largest_cluster_1.push_back(P30_largest_cluster_B / P30_largest_cluster_A.size());
        //std::cout << "P30_largest_cluster: " << P30_largest_cluster_B / nv << std::endl;

        double P32_largest_cluster_B = 0;
        for (size_t i = 0; i < P32_largest_cluster_A.size(); ++i)
        {
            P32_largest_cluster_B = P32_largest_cluster_B + P32_largest_cluster_A[i];
        }
        P32_largest_cluster_1.push_back(P32_largest_cluster_B / P32_largest_cluster_A.size());
        //std::cout << "P32_largest_cluster: " << P32_largest_cluster_B / nv << std::endl;

        double P30_connected_B = 0;
        for (size_t i = 0; i < P30_connected_A.size(); ++i)
        {
            P30_connected_B = P30_connected_B + P30_connected_A[i];
        }
        P30_connected_1.push_back(P30_connected_B / P30_connected_A.size());

        double Ratio_of_P30_B = 0;
        for (size_t i = 0; i < Ratio_of_P30_A.size(); ++i)
        {
            Ratio_of_P30_B = Ratio_of_P30_B + Ratio_of_P30_A[i];
        }
        Ratio_of_P30_1.push_back(Ratio_of_P30_B / Ratio_of_P30_A.size());

        for (size_t i = 0; i < Percolation_probability_A.size(); ++i)
        {
            if (Percolation_probability_A[i] == 1)
                nf++;
        }

        Percolation_probability_1.push_back(nf / ((double)(nv * 1.00)));
        double Percolation_probability_B = nf / ((double)(nv * 1.00));

        double n_I_B = 0;
        for (size_t i = 0; i < n_I_A.size(); ++i)
        {
            n_I_B = n_I_B + n_I_A[i];
        }
        n_I_1.push_back(n_I_B / n_I_A.size());

        if (str_frac_size == "powerlaw")
        {
            Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step.txt", array12[0][1], array12[0][2], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv, P30_connected_B / nv, Ratio_of_P30_B / nv);
        }
        else if (str_frac_size == "lognormal")
        {
            Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step.txt", array12[0][2], array12[0][3], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv, P30_connected_B / nv, Ratio_of_P30_B / nv);
        }
        else if (str_frac_size == "uniform")
        {
            Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step.txt", array12[0][0], array12[0][1], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv, P30_connected_B / nv, Ratio_of_P30_B / nv);
        }
        else if (str_frac_size == "single")
        {
            Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step.txt", array12[0][0], array12[0][0], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv, P30_connected_B / nv, Ratio_of_P30_B / nv);
        }

        if (njk == 0 && Percolation_probability_B > 0.49999)
        {
            njk++;
            Percolation_parameter_c = Percolation_parameter_B_1 / Percolation_parameter_A_1.size();
            std::cout << "\n**********Found Pc**********\n\n";
            //Sign_of_finding_pc("Pc_Found.txt");
        }

        if (njk != 0 && Percolation_parameter_B_1 / Percolation_parameter_A_1.size() > 2 * Percolation_parameter_c)
        {
            std::cout << "\n**********Found two times Pc**********\n\n";
            break;
        }
    };

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

inline void Loop_DFN::Sign_of_finding_pc(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "Pc has been found: ";
    oss << Percolation_parameter_c;

    oss.close();
}
}; // namespace DFN
