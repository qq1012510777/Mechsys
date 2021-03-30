#ifndef LOOP_H
#define LOOP_H

#include <mechsys/dfn_c/Domain.h>
#include <omp.h>

namespace DFN
{
    class Loop
    {
    public:
        double times;                   ///< loop times, each time DFN modeling, the density will be increased compared to last time DFN modeling
                                        //double RatioLR;                 ///< ratio of L and Ra, L is model edge size, Ra is the average radius of min and max ourter circle radius
                                        //double R_a;                     ///< the average radius of min and max ourter circle radius
                                        //double R_low;                   ///< minimum radius of outer circle radius
                                        //double R_up;                    ///< maximum radius of outer circle radius
        Array<Vec4_t> array12;          ///< alpha (power law), min_radius, max_radius,
        double L;                       ///< model size
        size_t Nproc;                   ///< num of threads
        size_t nt;                      ///< should not larger than times! which model is the model shown to us
        size_t nk;                      ///< when nk DFN models are finished, output one time
        size_t nv_MC_TIMES;             ///< each density, the MC times
        double nx;                      ///< the increment of fracture number regard to each DFN modeling time
        double Percolation_parameter_c; ///< when percolation probability equals to 0.5, the percolation parameter is
        size_t NumofFsets;              ///< number of fracture sets
        double P_end;                   ///< end condition of adding fractures
        double L_increment;             ///< model edge increment
        Array<double> DenWeight;        ///< the weight of number of each sets
        Array<Vec7_t> array13;          ///< the size of this array is the number of fracture sets; and each element includes the seven inputs of each Fisher fracture set, inputs are mean dip direction, mean dip angle, Fisher constant, min dip direction, max dip direction, min dip angle, max dip angle
        Array<double> P32_total_1;      ///< the next several 1D arrays store outputs
        Array<double> P32_connected_1;
        Array<double> P30_1;
        Array<double> Percolation_parameter_1;
        Array<double> Percolation_parameter_2;
        Array<double> Percolation_parameter_3;
        Array<double> Percolation_parameter_4;
        Array<double> Percolation_parameter_5;
        Array<double> Ratio_of_P32_1;
        Array<double> Percolation_probability_1;
        Array<double> n_I_1;
        Array<double> Correlation_length_1;
        Array<double> Max_gyration_radius_1;
        Array<double> P30_largest_cluster_1;
        Array<double> P32_largest_cluster_1;

    public:
        void Loop_create_DFNs(gsl_rng *random_seed, String str_ori, String str_frac_size, String percolation_direction);
        void Data_output_stepBYstep(size_t times, char const *FileKey, double min_R, double max_R, double L, double increment_fracture, double P32_total_B, double P32_connected_B, double P30_B, double Ratio_of_P32_B, double Percolation_parameter_B_1, double Percolation_parameter_B_2, double Percolation_parameter_B_3, double Percolation_parameter_B_4, double Percolation_parameter_B_5, double n_I, double Percolation_probability_B, double Correlation_length_B, double Gyration_radius_B, double P30_largest_cluster_B, double P32_largest_cluster_B);
        void Sign_of_finding_pc(char const *FileKey); ///< if Pc is found, outputs a file
    };

    void Loop::Loop_create_DFNs(gsl_rng *random_seed, String str_ori, String str_frac_size, String percolation_direction)
    {

        size_t nv = nv_MC_TIMES; //each density, the MC times

        double R_up;
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

        size_t np = 0;
        size_t njk = 0;

        while (np < times)
        {
            np++;
            size_t n = /*np **/ nx;

            double nf = 0;
            L = L + L_increment;

            Array<Vec2_t> array11; //min x, max x; min y, max y; min z, max z; (they are DOMAIN size, not MODEL size!!!)
            array11.resize(3);
            array11[0][0] = -L * 0.5;                                                         // - R_up;
            array11[0][1] = L * 0.5;                                                          // + R_up;
            array11[1][0] = -L * 0.5;                                                         // - R_up;
            array11[1][1] = L * 0.5;                                                          // + R_up;
            array11[2][0] = -L * 0.5;                                                         // - R_up;
            array11[2][1] = L * 0.5;                                                          // + R_up;
            double model_size[6] = {-L * 0.5, L * 0.5, -L * 0.5, L * 0.5, -L * 0.5, L * 0.5}; // MODEL size, not DOMAIN

            Array<double> P32_total_A;
            Array<double> P32_connected_A;
            Array<double> P30_A;
            Array<double> Percolation_parameter_A_1;
            Array<double> Percolation_parameter_A_2;
            Array<double> Percolation_parameter_A_3;
            Array<double> Percolation_parameter_A_4;
            Array<double> Percolation_parameter_A_5;
            Array<double> Ratio_of_P32_A;
            Array<double> Percolation_probability_A;
            Array<double> n_I_A;
            Array<double> Correlation_length_A;
            Array<double> Max_gyration_radius_A;
            Array<double> P30_largest_cluster_A;
            Array<double> P32_largest_cluster_A;

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

            size_t k_index = 0;

            DFN::Domain dom_1;
            dom_1.Create_whole_model(n, DenWeight, random_seed, model_size, str_ori, str_frac_size, array11, array12, array13, P_end,1);
            double P32_total_B;
            double P32_connected_B;
            double P30_B;
            double Percolation_parameter_B_1;
            double Percolation_parameter_B_2;
            double Percolation_parameter_B_3;
            double Percolation_parameter_B_4;
            double Percolation_parameter_B_5;
            double Ratio_of_P32_B;
            double Correlation_length_B;
            double Max_gyration_radius_B;
            double P30_largest_cluster_B;
            double P32_largest_cluster_B;
            double Percolation_probability_B;
            double n_I_B;

            if(dom_1.Percolation_parameter_a == -100)
            {
                k_index = -100;
                goto AH100;    
            }else
                n = dom_1.Fractures.Size();


#pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i = 0; i < nv; i++)
            {
                DFN::Domain dom;
                //std::cout<<"debug\n";
                size_t z;
                dom.Create_whole_model(n, DenWeight, random_seed, model_size, str_ori, str_frac_size, array11, array12, array13, P_end,0); ///uniform means oritation data are generated uniformly, so, actually, array13 is input but not used
                if (dom.Percolation_parameter_a == -100)
                {
                    k_index = -100;
                    goto k100;
                }
                

                z = dom.Identify_percolation_clusters(percolation_direction);
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

                if (z == 1)
                    Percolation_probability_A[i] = 1;
                else
                    Percolation_probability_A[i] = 0;

                if (np == nt && i == nv - 1)
                {
                    //dom.WriteFrac("tdfn01");
                    dom.PlotMatlab_DFN("tdfn01_DFN");
                    dom.PlotMatlab_DFN_and_Intersection("tdfn01_DFN_and_Intersections");
                    dom.PlotMatlab_ORI_SCATTER("tdfn01_ORI_SCATTER");
                    //dom.PlotMatlab_Traces_on_Model_surfaces("tdfn01_Trace_on_surfaces");
                    dom.PlotMatlab_DFN_Highlight_Cluster("tdfn01_DFN_Highlight_Cluster");
                    dom.PLotMatlab_DFN_Cluster_along_a_direction("tdfn01_DFN_Z_clusters", "z");
                    dom.PlotMatlab_Radius_and_Area_kstest("tdfn01_DFN_Fracture_Radius_and_Area");
                    dom.PlotMatlab_Radius_and_Perimeter_kstest("tdfn01_DFN_Fracture_Radius_and_Perimeter");
                    dom.DataFile_Radius_AreaAndPerimeter("tdfn01_DFN_Radius_AreaAndPerimeter");
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
            k100:;
            }

            if (k_index == -100)
                goto AH100;


            P32_total_B = 0;
            for (size_t i = 0; i < P32_total_A.size(); ++i)
            {
                P32_total_B = P32_total_B + P32_total_A[i];
            }
            P32_total_1.Push(P32_total_B / P32_total_A.size());

            P32_connected_B = 0;
            for (size_t i = 0; i < P32_connected_A.size(); ++i)
            {
                P32_connected_B = P32_connected_B + P32_connected_A[i];
            }
            P32_connected_1.Push(P32_connected_B / P32_connected_A.size());

            P30_B = 0;
            for (size_t i = 0; i < P30_A.size(); ++i)
            {
                P30_B = P30_B + P30_A[i];
            }
            P30_1.Push(P30_B / P30_A.size());

            Percolation_parameter_B_1 = 0;
            for (size_t i = 0; i < Percolation_parameter_A_1.size(); ++i)
            {
                Percolation_parameter_B_1 = Percolation_parameter_B_1 + Percolation_parameter_A_1[i];
            }
            Percolation_parameter_1.Push(Percolation_parameter_B_1 / Percolation_parameter_A_1.size());

            Percolation_parameter_B_2 = 0;
            for (size_t i = 0; i < Percolation_parameter_A_2.size(); ++i)
            {
                Percolation_parameter_B_2 = Percolation_parameter_B_2 + Percolation_parameter_A_2[i];
            }
            Percolation_parameter_2.Push(Percolation_parameter_B_2 / Percolation_parameter_A_2.size());

            Percolation_parameter_B_3 = 0;
            for (size_t i = 0; i < Percolation_parameter_A_3.size(); ++i)
            {
                Percolation_parameter_B_3 = Percolation_parameter_B_3 + Percolation_parameter_A_3[i];
            }
            Percolation_parameter_3.Push(Percolation_parameter_B_3 / Percolation_parameter_A_3.size());

            Percolation_parameter_B_4 = 0;
            for (size_t i = 0; i < Percolation_parameter_A_4.size(); ++i)
            {
                Percolation_parameter_B_4 = Percolation_parameter_B_4 + Percolation_parameter_A_4[i];
            }
            Percolation_parameter_4.Push(Percolation_parameter_B_4 / Percolation_parameter_A_4.size());

            Percolation_parameter_B_5 = 0;
            for (size_t i = 0; i < Percolation_parameter_A_5.size(); ++i)
            {
                Percolation_parameter_B_5 = Percolation_parameter_B_5 + Percolation_parameter_A_5[i];
            }
            Percolation_parameter_5.Push(Percolation_parameter_B_5 / Percolation_parameter_A_5.size());

            Ratio_of_P32_B = 0;
            for (size_t i = 0; i < Ratio_of_P32_A.size(); ++i)
            {
                Ratio_of_P32_B = Ratio_of_P32_B + Ratio_of_P32_A[i];
            }
            Ratio_of_P32_1.Push(Ratio_of_P32_B / Ratio_of_P32_A.size());

            Correlation_length_B = 0;
            for (size_t i = 0; i < Correlation_length_A.size(); ++i)
            {
                Correlation_length_B = Correlation_length_B + Correlation_length_A[i];
            }
            Correlation_length_1.Push(Correlation_length_B / Correlation_length_A.size());

            Max_gyration_radius_B = 0;
            for (size_t i = 0; i < Max_gyration_radius_A.size(); ++i)
            {
                Max_gyration_radius_B = Max_gyration_radius_B + Max_gyration_radius_A[i];
            }
            Max_gyration_radius_1.Push(Max_gyration_radius_B / Max_gyration_radius_A.size());

            P30_largest_cluster_B = 0;
            for (size_t i = 0; i < P30_largest_cluster_A.size(); ++i)
            {
                P30_largest_cluster_B = P30_largest_cluster_B + P30_largest_cluster_A[i];
            }
            P30_largest_cluster_1.Push(P30_largest_cluster_B / P30_largest_cluster_A.size());
            //std::cout << "P30_largest_cluster: " << P30_largest_cluster_B / nv << std::endl;

            P32_largest_cluster_B = 0;
            for (size_t i = 0; i < P32_largest_cluster_A.size(); ++i)
            {
                P32_largest_cluster_B = P32_largest_cluster_B + P32_largest_cluster_A[i];
            }
            P32_largest_cluster_1.Push(P32_largest_cluster_B / P32_largest_cluster_A.size());
            //std::cout << "P32_largest_cluster: " << P32_largest_cluster_B / nv << std::endl;

            for (size_t i = 0; i < Percolation_probability_A.size(); ++i)
            {
                if (Percolation_probability_A[i] == 1)
                    nf++;
            }
            Percolation_probability_1.Push(nf / ((double)(nv * 1.00)));
            Percolation_probability_B = nf / ((double)(nv * 1.00));

            n_I_B = 0;
            for (size_t i = 0; i < n_I_A.size(); ++i)
            {
                n_I_B = n_I_B + n_I_A[i];
            }
            n_I_1.Push(n_I_B / n_I_A.size());

        AH100:;
            if (k_index == -100 && np != 1)
                continue;
            if (str_frac_size == "powerlaw")
            {
                Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step", array12[0][1], array12[0][2], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv);
            }
            else if (str_frac_size == "lognormal")
            {
                Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step", array12[0][2], array12[0][3], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv);
            }
            else if (str_frac_size == "uniform")
            {
                Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step", array12[0][0], array12[0][1], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv);
            }
            else if (str_frac_size == "single")
            {
                Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step", array12[0][0], array12[0][0], L, nx, P32_total_B / nv, P32_connected_B / nv, P30_B / nv, Ratio_of_P32_B / nv, Percolation_parameter_B_1 / nv, Percolation_parameter_B_2 / nv, Percolation_parameter_B_3 / nv, Percolation_parameter_B_4 / nv, Percolation_parameter_B_5 / nv, n_I_B / n_I_A.size(), Percolation_probability_B, Correlation_length_B / nv, Max_gyration_radius_B / nv, P30_largest_cluster_B / nv, P32_largest_cluster_B / nv);
            }
            /*
            if (njk == 0 && Percolation_probability_B > 0.49999)
            {
                njk++;
                Percolation_parameter_c = Percolation_parameter_B_1 / Percolation_parameter_A_1.size();
                std::cout << "\n**********Found Pc**********\n\n";
                Sign_of_finding_pc("Pc_Found");
            }

            if (njk != 0 && Percolation_parameter_B_1 / Percolation_parameter_A_1.size() > 2 * Percolation_parameter_c)
            {
                std::cout << "\n**********Found two times Pc**********\n\n";
                break;
            }
            */
        }
        PLotMatlab_DFN_Connectivity("tdfn01_Connectivity", Percolation_parameter_1, Percolation_parameter_2, Percolation_parameter_3, Percolation_parameter_4, Percolation_parameter_5, n_I_1, P30_1, P32_total_1, Percolation_probability_1, array12[0][0]);
        PLotMatlab_DFN_Correlation_length_over_model_size("tdfn01_Correlation_length_over_model_size", Percolation_parameter_1, Percolation_parameter_2, Percolation_parameter_3, Percolation_parameter_4, Percolation_parameter_5, n_I_1, P30_1, P32_total_1, Correlation_length_1, L);
        PLotMatlab_DFN_max_gyration_radius_over_model_size("tdfn01_max_gyration_radius_over_model_size", Percolation_parameter_1, Percolation_parameter_2, Percolation_parameter_3, Percolation_parameter_4, Percolation_parameter_5, n_I_1, P30_1, P32_total_1, Max_gyration_radius_1, L);
        std::cout << "Loop finished!\n";
        //Datafile_output("tdfn01_Datafile", model_side_length, P32_total_1, P32_connected_1, P30_1, Percolation_parameter_1, Ratio_of_P32_1, Percolation_probability_1);
    }; // namespace DFN

    inline void Loop::Data_output_stepBYstep(size_t times, char const *FileKey, double min_R, double max_R, double L, double increment_fracture, double P32_total_B, double P32_connected_B, double P30_B, double Ratio_of_P32_B, double Percolation_parameter_B_1, double Percolation_parameter_B_2, double Percolation_parameter_B_3, double Percolation_parameter_B_4, double Percolation_parameter_B_5, double n_I, double Percolation_probability_B, double Correlation_length_B, double Gyration_radius_B, double P30_largest_cluster_B, double P32_largest_cluster_B)
    {
        if (times == 1)
        {
            std::ostringstream oss;
            oss << "Initial_model_edge"
                << "\t" << L << "\n";
            oss << "Fracture_increment:"
                << "\t" << increment_fracture << "\n";
            oss << "min_R"
                << "\t" << min_R << "\n";
            oss << "max_R"
                << "\t" << max_R << "\n";
            oss << "Model_edge"
                << "\t"
                << "P32_total"
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
                << "\n";
            oss << L << "\t" << P32_total_B << "\t" << P32_connected_B << "\t" << P30_B << "\t" << Ratio_of_P32_B << "\t" << Percolation_parameter_B_1 << "\t" << Percolation_parameter_B_2 << "\t" << Percolation_parameter_B_3 << "\t" << Percolation_parameter_B_4 << "\t" << Percolation_parameter_B_5 << "\t" << n_I << "\t" << Percolation_probability_B << "\t" << Correlation_length_B / L << "\t" << Gyration_radius_B / L << "\t" << P30_largest_cluster_B << "\t" << P32_largest_cluster_B << "\n";
            String fn(FileKey);
            fn.append(".txt");
            std::ofstream of(fn.CStr(), std::ios::out);
            of << oss.str();
            of.close();
        }
        else
        {
            std::ostringstream oss;
            oss << L << "\t" << P32_total_B << "\t" << P32_connected_B << "\t" << P30_B << "\t" << Ratio_of_P32_B << "\t" << Percolation_parameter_B_1 << "\t" << Percolation_parameter_B_2 << "\t" << Percolation_parameter_B_3 << "\t" << Percolation_parameter_B_4 << "\t" << Percolation_parameter_B_5 << "\t" << n_I << "\t" << Percolation_probability_B << "\t" << Correlation_length_B / L << "\t" << Gyration_radius_B / L << "\t" << P30_largest_cluster_B << "\t" << P32_largest_cluster_B << "\n";
            String fn(FileKey);
            fn.append(".txt");
            std::ofstream of(fn.CStr(), std::ios::app);
            of << oss.str();
            of.close();
        }
    };

    void Loop::Sign_of_finding_pc(char const *FileKey)
    {
        std::ostringstream oss;
        oss << "Pc has been found: ";
        oss << Percolation_parameter_c;

        String fn(FileKey);
        fn.append(".txt");
        std::ofstream of(fn.CStr(), std::ios::out);
        of << oss.str();
        of.close();
    }
}; // namespace DFN
#endif
