#pragma once
//#include "../FEM_H/FEM_DFN_A.h"
#include "../FEM_H/MHFEM.h"
#include "../HDF5_DFN/HDF5_DFN.h"
#include "../MATLAB_DATA_API/MATLAB_DATA_API.h"
#include "../Mesh_H/Mesh_DFN_linear.h"
#include "../ProgressBar/ProgressBar.h"
#include "Domain_WL.h"
#include "Load_a_DFN_from_matfile.h"
#include <cstdio>
#include <omp.h>
#include <stdio.h>
#include <unistd.h>

namespace DFN
{
class Multi_processes
{
public:
    double times;                  ///< loop times, each time DFN modeling, the density will be increased compared to last time DFN modeling
    std::vector<Vector4d> array12; ///< alpha (power law), min_radius, max_radius,
    double L;                      ///< model size
    size_t Nproc;                  ///< num of threads for DFN modeling
    size_t Nproc_flow;             ///< num of threads for flow modeling
    size_t nt;                     ///< should not larger than times! which model is the model shown to us
    size_t nk;                     ///< when nk DFN models are finished, output one time
    size_t nv_MC_TIMES;            ///< each density, the MC times
    double nx;                     ///< the increment of fracture number regard to each DFN modeling time
    double Density_c;              ///< when percolation probability equals to 0.5, the densit is
    size_t NumofFsets;             ///< number of fracture sets
    size_t Nb_flow_sim_MC_times;

    size_t n_initial_frac_density = 0;

    std::vector<double> DenWeight; ///< the weight of number of each sets
    std::vector<Vector7d> array13; ///< the size of this array is the number of fracture sets;
    //and each element includes the seven inputs of each Fisher fracture set,
    //inputs are mean dip direction, mean dip angle, Fisher constant, min dip direction, max dip direction, min dip angle, max dip angle

    string Data_CommandFile = "Data_Command.m";
    string Data_MatFile = "DFN_file.h5";

    size_t switch_2D = 0;

    size_t Model_flow = 0;

public:
    void Loop_create_DFNs(gsl_rng *random_seed,
                          string str_ori,
                          string str_frac_size,
                          string percolation_direction,
                          const double min_ele_edge,
                          const double max_ele_edge,
                          const string conductivity_distri,
                          size_t modelno);

    void Sign_of_finding_pc(string FileKey);
    ///< if Pc is found, outputs a file

    void Matlab_Data_output_stepBYstep(const size_t np,
                                       string FileKey_mat,
                                       std::vector<double> P32_total_A,
                                       std::vector<std::vector<double>> P32_connected_A,
                                       std::vector<double> P30_A,
                                       std::vector<std::vector<double>> Ratio_of_P32_A,
                                       std::vector<std::vector<double>> Percolation_probability_A,
                                       std::vector<double> n_I_A,
                                       std::vector<double> P30_largest_cluster_A,
                                       std::vector<double> P32_largest_cluster_A,
                                       std::vector<std::vector<double>> P30_connected_A,
                                       std::vector<std::vector<double>> Ratio_of_P30_A,
                                       std::vector<std::vector<double>> Permeability_A,
                                       std::vector<std::vector<double>> Q_error_A,
                                       const string str_ori,
                                       const string str_frac_size,
                                       const string conductivity_distri,
                                       const double domain_size,
                                       std::vector<std::vector<double>> min_ele_edge,
                                       std::vector<std::vector<double>> max_ele_edge,
                                       std::vector<std::vector<double>> Mean_ele_edge);

    void Matlab_command(string FileKey_m, string FileKey_mat, size_t np, size_t ny, size_t model_no);

private:
    Vector2d Re_MHFEM(DFN::Domain dom_222,
                      const double min_ele_edge,
                      const double max_ele_edge,
                      size_t dir_,
                      double h_in,
                      double h_out,
                      size_t Nproc_t_);
};

//****************************
inline void Multi_processes::Loop_create_DFNs(gsl_rng *random_seed,
                                              string str_ori,
                                              string str_frac_size,
                                              string percolation_direction,
                                              const double min_ele_edge,
                                              const double max_ele_edge,
                                              const string conductivity_distri,
                                              size_t modelno)
{
    if (this->Model_flow == 1)
    {
        cout << "\033[31m------ switch the flow modeling on ------";
        cout << "\033[0m" << endl;
    }
    //bool if_probability_1 = false;
    //bool if_probability_2 = false;
    //bool show_flow = false;

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

    //string DFN_h5_file = "DFN_geom_" + To_string_with_width(modelno, 3) + ".h5";
    //const char *CS = DFN_h5_file.c_str();
    //std::remove(CS);
    //DFN::HDF5_DFN DFN_geome(DFN_h5_file);

    while (np < times)
    {
        np++;
        size_t n = np * nx;

        //double nf = 0;

        std::vector<double> P32_total_A;
        std::vector<std::vector<double>> P32_connected_A(3); // x, y, z
        std::vector<double> P30_A;
        std::vector<std::vector<double>> Ratio_of_P32_A(3);            // x, y, z
        std::vector<std::vector<double>> Percolation_probability_A(3); // x, y, z
        std::vector<double> n_I_A;
        std::vector<double> P30_largest_cluster_A;
        std::vector<double> P32_largest_cluster_A;
        std::vector<std::vector<double>> P30_connected_A(3);
        std::vector<std::vector<double>> Ratio_of_P30_A(3);
        std::vector<std::vector<double>> Permeability_A(3);
        std::vector<std::vector<double>> Q_error_A(3);
        std::vector<DFN::Domain> Dom_vec;
        std::vector<std::vector<DFN::Mesh_DFN_linear>> Mesh_vec(3);
        std::vector<std::vector<size_t>> Mesh_status(3);

        std::vector<std::vector<double>> min_edge_r(3);
        std::vector<std::vector<double>> max_edge_r(3);
        std::vector<std::vector<double>> mean_edge_r(3);

        // 0: no percolating cluster
        // 1: have percolation cluster, and mesh successfully
        // 2: have percolation cluster, but fail to mesh

        P32_total_A.resize(nv);
        P32_connected_A[0].resize(nv);
        P32_connected_A[1].resize(nv);
        P32_connected_A[2].resize(nv);
        P30_A.resize(nv);
        Ratio_of_P32_A[0].resize(nv);
        Ratio_of_P32_A[1].resize(nv);
        Ratio_of_P32_A[2].resize(nv);
        Percolation_probability_A[0].resize(nv);
        Percolation_probability_A[1].resize(nv);
        Percolation_probability_A[2].resize(nv);
        n_I_A.resize(nv);
        P30_largest_cluster_A.resize(nv);
        P32_largest_cluster_A.resize(nv);
        P30_connected_A[0].resize(nv);
        P30_connected_A[1].resize(nv);
        P30_connected_A[2].resize(nv);
        Ratio_of_P30_A[0].resize(nv);
        Ratio_of_P30_A[1].resize(nv);
        Ratio_of_P30_A[2].resize(nv);

        min_edge_r[0].resize(nv);
        min_edge_r[1].resize(nv);
        min_edge_r[2].resize(nv);

        max_edge_r[0].resize(nv);
        max_edge_r[1].resize(nv);
        max_edge_r[2].resize(nv);

        mean_edge_r[0].resize(nv);
        mean_edge_r[1].resize(nv);
        mean_edge_r[2].resize(nv);

        if (nv < Nb_flow_sim_MC_times)
            Nb_flow_sim_MC_times = nv;

        //--------------
        Dom_vec.resize(Nb_flow_sim_MC_times);
        Mesh_vec[0].resize(Nb_flow_sim_MC_times);
        Mesh_vec[1].resize(Nb_flow_sim_MC_times);
        Mesh_vec[2].resize(Nb_flow_sim_MC_times);

        Mesh_status[0] = vector<size_t>(Nb_flow_sim_MC_times, 0);
        Mesh_status[1] = vector<size_t>(Nb_flow_sim_MC_times, 0);
        Mesh_status[2] = vector<size_t>(Nb_flow_sim_MC_times, 0);

        Permeability_A[0].resize(Nb_flow_sim_MC_times);
        Permeability_A[1].resize(Nb_flow_sim_MC_times);
        Permeability_A[2].resize(Nb_flow_sim_MC_times);

        Q_error_A[0].resize(Nb_flow_sim_MC_times);
        Q_error_A[1].resize(Nb_flow_sim_MC_times);
        Q_error_A[2].resize(Nb_flow_sim_MC_times);

        //DFN::ProgressBar prog_bar;
        std::cout << "\nThe Model NO." << np << " is creating now! Sizes: " << str_frac_size << "; Ori: " << str_ori << ".\n";
    Regenerate_dfn:;
        nv = 1;
        for (size_t i = 0; i < nv; i++)
        {
            try
            {
                DFN::Domain dom;
                //dom.Fractures.clear();
                if (switch_2D == 0)
                {
                    dom.Create_whole_model((n + n_initial_frac_density),
                                           DenWeight,
                                           random_seed,
                                           model_size,
                                           str_ori,
                                           str_frac_size,
                                           array11,
                                           array12,
                                           array13,
                                           conductivity_distri);
                }
                else if (switch_2D == 1)
                {

                    dom.mode_2D = true;

                    dom.Create_whole_model((n + n_initial_frac_density),
                                           DenWeight,
                                           random_seed,
                                           model_size,
                                           str_ori,
                                           str_frac_size,
                                           array11,
                                           array12,
                                           array13,
                                           conductivity_distri);
                }
                else
                {
                    throw Error_throw_pause("Undefined mode!\n");
                }
                ///uniform means oritation data
                //are generated uniformly, so,
                //actually, array13 is input but not used
                if (switch_2D == 1 && percolation_direction != "z")
                    throw Error_throw_pause("When it is in 2D mode, percolation direction must be along z direction!\n");

                dom.Identify_percolation_clusters();

                dom.Connectivity_analysis();

                P32_total_A[i] = (dom.P32_total);
                P32_connected_A[0][i] = (dom.P32_connected[0]);
                P32_connected_A[1][i] = (dom.P32_connected[1]);
                P32_connected_A[2][i] = (dom.P32_connected[2]);
                P30_A[i] = (dom.P30);
                Ratio_of_P32_A[0][i] = (dom.Ratio_of_P32[0]);
                Ratio_of_P32_A[1][i] = (dom.Ratio_of_P32[1]);
                Ratio_of_P32_A[2][i] = (dom.Ratio_of_P32[2]);
                n_I_A[i] = (dom.n_I);
                P30_largest_cluster_A[i] = (dom.P30_largest_cluster);
                P32_largest_cluster_A[i] = (dom.P32_largest_cluster);
                P30_connected_A[0][i] = (dom.P30_connected[0]);
                P30_connected_A[1][i] = (dom.P30_connected[1]);
                P30_connected_A[2][i] = (dom.P30_connected[2]);
                Ratio_of_P30_A[0][i] = (dom.Ratio_of_P30[0]);
                Ratio_of_P30_A[1][i] = (dom.Ratio_of_P30[1]);
                Ratio_of_P30_A[2][i] = (dom.Ratio_of_P30[2]);

                if (dom.Percolation_status[0] == true)
                    Percolation_probability_A[0][i] = 1;
                else
                    Percolation_probability_A[0][i] = 0;

                if (dom.Percolation_status[1] == true)
                    Percolation_probability_A[1][i] = 1;
                else
                    Percolation_probability_A[1][i] = 0;

                if (dom.Percolation_status[2] == true)
                    Percolation_probability_A[2][i] = 1;
                else
                    Percolation_probability_A[2][i] = 0;

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
                    dom.PlotMatlab_Radius_and_Perimeter("tdfn01_DFN_Fracture_Radius_and_Perimeter.m");
                    //dom.DataFile_Radius_AreaAndPerimeter("tdfn01_DFN_Radius_AreaAndPerimeter.m");
                    dom.Matlab_Out_Frac_matfile("Fractures.mat");
                }

                if (i < Nb_flow_sim_MC_times /*&& this->Model_flow == 1*/)
                {
                    Dom_vec[i] = dom;
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
                //exit(0);
            }
            catch (bad_alloc &e)
            {
                cout << "\033[33mRegenerate a DFN! Because:\n"
                     << e.what() << "\033[0m" << endl;
                goto Regenerate_dfn;
            }
            catch (...)
            {
                goto Regenerate_dfn;
            }

            //prog_bar.Rep_prog_for_paral(nv, 5, "\t\tDFN_MC_modeling ");
        }
        //cout << endl;
        // for-loop for connectivity ends here

        /*
        for (size_t i = 0; i < Nb_flow_sim_MC_times; ++i)
        {
            DFN::Domain dom = Dom_vec[i];
            string groupname_ = "group_" + To_string_with_width(np, 2) + "_MCtime_" + To_string_with_width(i + 1, 5);

            dom.Output_all_fractures_HDF5(DFN_h5_file, groupname_, 5);

            if (i == 0)
                DFN_geome.Append_dataset_to_group(DFN_h5_file, groupname_, "MC_times", vector<double>{(double)Nb_flow_sim_MC_times});
        }
        */

        if (this->Model_flow == 1)
        {
            //DFN::ProgressBar prog_bar_2;
            auto start_1 = std::chrono::steady_clock::now();
            cout << "\033[33m\tmeshing started!\033[0m" << endl;
            for (size_t i = 0; i < Nb_flow_sim_MC_times; ++i)
            {
                DFN::Domain dom;
                dom = Dom_vec[i];

                for (size_t rt = 2; rt < 3; ++rt)
                {
                    try
                    {
                        //cout << "flow model " << rt << endl;
                        bool z = dom.Percolation_status[rt];
                        if (z == true)
                        {
                            dom.Re_identify_intersection_considering_trimmed_frac();
                            dom.Identify_percolation_clusters();
                            //dom.PlotMatlab_DFN_and_Intersection("tdfn01_DFN_and_Intersections_II.m");
                            if (dom.Percolation_status[rt] != z)
                            {
                                throw Error_throw_ignore("Trimmed fractures did not form at least a percolation cluster!\n");
                            }

                            //cout << "\t\tmesh start" << endl;
                            DFN::Mesh_DFN_linear mesh_s(dom, min_ele_edge, max_ele_edge, rt, Nproc);
                            //cout << "\t\tmesh finish" << endl;

                            Mesh_vec[rt][i] = mesh_s;
                            Mesh_status[rt][i] = 1;

                            min_edge_r[rt][i] = mesh_s.min_edge_length;
                            max_edge_r[rt][i] = mesh_s.max_edge_length;
                            mean_edge_r[rt][i] = mesh_s.mean_edge_length;

                            if (mesh_s.mesh_state == false)
                            {
                                Mesh_status[rt][i] = 2;
                                goto Regenerate_dfn;
                            }
                        }
                        else
                        {
                            Mesh_status[rt][i] = 0;
                        }
                    }
                    catch (std::string const &ErrMsg)
                    {
                        //cout << "\t\tmesh error! label this mesh as failed one" << endl;
                        Mesh_status[rt][i] = 2;
                        goto Regenerate_dfn;
                    }
                    catch (...)
                    {
                        //cout << "\t\tOther error! label this mesh as failed one" << endl;
                        Mesh_status[rt][i] = 2;
                        goto Regenerate_dfn;
                    }
                }

                //_2.Rep_prog_serially_for_supercomputer(i + 1, 5, Nb_flow_sim_MC_times, "\t\tMeshing ");
            }
            //cout << endl;
            auto end_1 = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::micro> elapsed_1 = end_1 - start_1; // std::micro time (us)
            cout << "\033[33m\tmeshing finished! runtime: " << (((double)(elapsed_1.count() * 1.0) * (0.000001)) / 60.00) << " minutes \033[0m" << endl;
        };

        if (this->Model_flow == 1)
        {
            //DFN::ProgressBar prog_bar_2;

            auto start_2 = std::chrono::steady_clock::now();
            cout << "\033[31m\tFEM started!\033[0m" << endl;

            for (size_t i = 0; i < Nb_flow_sim_MC_times; ++i)
            {
                DFN::Domain dom;
                dom = Dom_vec[i];

                for (size_t rt = 2; rt < 3; ++rt)
                {
                    try
                    {
                        if (Mesh_status[rt][i] == 1)
                        {
                            DFN::Mesh_DFN_linear mesh_s = Mesh_vec[rt][i];

                            double h_in = 100, h_out = 20;
                            //cout << "\t\tFEM start" << endl;
                            DFN::MHFEM CC(mesh_s, dom, h_in, h_out, rt, this->Nproc_flow);
                            //cout << "\t\tFEM finish" << endl;
                            if (np == nt)
                            {
                                string mesh_ = "mesh_DFN_" + to_string(rt);
                                string FEM_ = "FEM_DFN_" + to_string(rt);
                                mesh_s.Matlab_plot(mesh_ + ".mat", mesh_ + ".m", dom);
                                CC.Matlab_plot(FEM_ + ".mat", FEM_ + ".m", mesh_s, dom);
                            }
                            Permeability_A[rt][i] = CC.Permeability;
                            Q_error_A[rt][i] = CC.Q_error;
                            //Permeability_A[rt][i] = 0;

                            if (CC.Q_error > 1 || isnan(CC.Permeability) == 1)
                            {
                                //cout << "\t\tfound large error or isnan, the error: " << CC.Q_error << ", the permeability: " << CC.Permeability << endl;
                                //dom.Matlab_Out_Frac_matfile("Fractures.mat");
                                //CC.Matlab_plot("MHFEM.mat", "MHFEM.m", mesh_s, dom);
                                //exit(0);
                                Vector2d TTE = this->Re_MHFEM(dom, min_ele_edge, max_ele_edge, rt, h_in, h_out, this->Nproc_flow);
                                Permeability_A[rt][i] = TTE[0];
                                Q_error_A[rt][i] = TTE[1];
                                if (TTE[0] > 1 || isnan(TTE[1]) == 1)
                                    goto Regenerate_dfn;
                            }
                        }
                        else if (Mesh_status[rt][i] == 2)
                        {
                            Permeability_A[rt][i] = -1;
                            Q_error_A[rt][i] = -1;
                        }
                        else if (Mesh_status[rt][i] == 0)
                        {
                            Permeability_A[rt][i] = 0;
                            Q_error_A[rt][i] = 0;
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
                        Permeability_A[rt][i] = -1;
                        Q_error_A[rt][i] = -1;
                        goto Regenerate_dfn;
                    }
                    catch (bad_alloc &e)
                    {
                        Permeability_A[rt][i] = -1;
                        Q_error_A[rt][i] = -1;
                        cout << "\t\tbad alloc\n";
                        goto Regenerate_dfn;
                    }
                    catch (...)
                    {
                        goto Regenerate_dfn;
                    }
                }

                //_2.Rep_prog_for_paral(Nb_flow_sim_MC_times, 5, "\t\tFEM ");
            }
            //cout << endl;
            auto end_2 = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::micro> elapsed_2 = end_2 - start_2; // std::micro time (us)
            cout << "\033[31m\tFEM finished! runtime: " << (((double)(elapsed_2.count() * 1.0) * (0.000001)) / 60.00) << " minutes \033[0m" << endl;
        }

        //---------------------------
        //cout << "\tstart output data\n";
        this->Matlab_Data_output_stepBYstep(np,
                                            Data_MatFile,
                                            P32_total_A,
                                            P32_connected_A,
                                            P30_A,
                                            Ratio_of_P32_A,
                                            Percolation_probability_A,
                                            n_I_A,
                                            P30_largest_cluster_A,
                                            P32_largest_cluster_A,
                                            P30_connected_A,
                                            Ratio_of_P30_A,
                                            Permeability_A,
                                            Q_error_A,
                                            str_ori,
                                            str_frac_size,
                                            conductivity_distri,
                                            this->L,
                                            min_edge_r,
                                            max_edge_r,
                                            mean_edge_r);
        //cout << "\tfinish output data\n\n";
    };

    //this->Matlab_command(Data_CommandFile, Data_MatFile, np, np, modelno);
    std::cout << "Loop finished!\n";
};

inline void Multi_processes::Matlab_Data_output_stepBYstep(const size_t np,
                                                           string FileKey_mat,
                                                           std::vector<double> P32_total_A,
                                                           std::vector<std::vector<double>> P32_connected_A,
                                                           std::vector<double> P30_A,
                                                           std::vector<std::vector<double>> Ratio_of_P32_A,
                                                           std::vector<std::vector<double>> Percolation_probability_A,
                                                           std::vector<double> n_I_A,
                                                           std::vector<double> P30_largest_cluster_A,
                                                           std::vector<double> P32_largest_cluster_A,
                                                           std::vector<std::vector<double>> P30_connected_A,
                                                           std::vector<std::vector<double>> Ratio_of_P30_A,
                                                           std::vector<std::vector<double>> Permeability_A,
                                                           std::vector<std::vector<double>> Q_error_A,
                                                           const string str_ori,
                                                           const string str_frac_size,
                                                           const string conductivity_distri,
                                                           const double domain_size,
                                                           vector<std::vector<double>> min_ele_edge,
                                                           vector<std::vector<double>> max_ele_edge,
                                                           vector<std::vector<double>> Mean_ele_edge)

{

    if (np == 1)
    {
        const char *CS = FileKey_mat.c_str();
        std::remove(CS);
        DFN::HDF5_DFN H5file_1(FileKey_mat); //create the file;

        H5file_1.Write_H5(FileKey_mat, "Ori_distr", str_ori);
        H5file_1.Write_H5(FileKey_mat, "FracSize_distr", str_frac_size);
        H5file_1.Write_H5(FileKey_mat, "Conduc_distr", conductivity_distri);
        H5file_1.Write_H5(FileKey_mat, "Domain_size", domain_size);

        H5file_1.Write_H5(FileKey_mat, "Frac_sets", (double)this->array12.size());

        for (size_t i = 0; i < this->array12.size(); ++i)
        {
            //cout << this->array12[i].rows() << endl;
            vector<double> tmp_1(this->array12[i].rows());

            for (size_t j = 0; j < (size_t)this->array12[i].rows(); ++j)
                tmp_1[j] = this->array12[i][j];

            string datasetname = "Frac_set_" + To_string_with_width(i + 1, 3) + "_sizes";

            if (i == 0)
                H5file_1.Write_H5(FileKey_mat, "Frac_sizes_oris", datasetname, tmp_1);
            else
                H5file_1.Append_dataset_to_group(FileKey_mat, "Frac_sizes_oris", datasetname, tmp_1);
        }

        if (str_ori == "fisher")
        {
            for (size_t i = 0; i < this->array13.size(); ++i)
            {
                //cout << this->array12[i].rows() << endl;
                vector<double> tmp_1(this->array13[i].cols());

                for (size_t j = 0; j < (size_t)this->array13[i].cols(); ++j)
                    tmp_1[j] = this->array13[i][j];

                string datasetname = "Frac_set_" + To_string_with_width(i + 1, 3) + "_oris";

                H5file_1.Append_dataset_to_group(FileKey_mat, "Frac_sizes_oris", datasetname, tmp_1);
            }
        }
    }
    DFN::HDF5_DFN H5file;

    string groupname = "group_" + To_string_with_width(np, 3);

    H5file.Overwrite(FileKey_mat, "Loop_times", np);

    string string_P32_total = "P32_total";
    string string_P32_connected_x = "P32_connected_x";
    string string_P32_connected_y = "P32_connected_y";
    string string_P32_connected_z = "P32_connected_z";
    string string_P30 = "P30";
    string string_Ratio_of_P32_x = "Ratio_of_P32_x";
    string string_Ratio_of_P32_y = "Ratio_of_P32_y";
    string string_Ratio_of_P32_z = "Ratio_of_P32_z";
    string string_Percolation_probability_x = "Percolation_probability_x";
    string string_Percolation_probability_y = "Percolation_probability_y";
    string string_Percolation_probability_z = "Percolation_probability_z";
    string string_n_I = "n_I";
    string string_P30_largest_cluster = "P30_largest_cluster";
    string string_P32_largest_cluster = "P32_largest_cluster";
    string string_P30_connected_x = "P30_connected_x";
    string string_P30_connected_y = "P30_connected_y";
    string string_P30_connected_z = "P30_connected_z";
    string string_Ratio_of_P30_x = "Ratio_of_P30_x";
    string string_Ratio_of_P30_y = "Ratio_of_P30_y";
    string string_Ratio_of_P30_z = "Ratio_of_P30_z";
    string string_Permeability_x = "Permeability_x";
    string string_Permeability_y = "Permeability_y";
    string string_Permeability_z = "Permeability_z";
    string string_Q_error_x = "Q_error_x";
    string string_Q_error_y = "Q_error_y";
    string string_Q_error_z = "Q_error_z";
    string string_Min_ele_edge_x = "Min_ele_edge_x";
    string string_Max_ele_edge_x = "Max_ele_edge_x";
    string string_Mean_ele_edge_x = "Mean_ele_edge_x";
    string string_Min_ele_edge_y = "Min_ele_edge_y";
    string string_Max_ele_edge_y = "Max_ele_edge_y";
    string string_Mean_ele_edge_y = "Mean_ele_edge_y";
    string string_Min_ele_edge_z = "Min_ele_edge_z";
    string string_Max_ele_edge_z = "Max_ele_edge_z";
    string string_Mean_ele_edge_z = "Mean_ele_edge_z";

    H5file.Write_H5(FileKey_mat, groupname, string_P32_total, P32_total_A); // created a group already
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P32_connected_x, P32_connected_A[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P32_connected_y, P32_connected_A[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P32_connected_z, P32_connected_A[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P30, P30_A);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Ratio_of_P32_x, Ratio_of_P32_A[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Ratio_of_P32_y, Ratio_of_P32_A[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Ratio_of_P32_z, Ratio_of_P32_A[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Percolation_probability_x, Percolation_probability_A[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Percolation_probability_y, Percolation_probability_A[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Percolation_probability_z, Percolation_probability_A[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_n_I, n_I_A);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P30_largest_cluster, P30_largest_cluster_A);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P32_largest_cluster, P32_largest_cluster_A);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P30_connected_x, P30_connected_A[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P30_connected_y, P30_connected_A[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_P30_connected_z, P30_connected_A[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Ratio_of_P30_x, Ratio_of_P30_A[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Ratio_of_P30_y, Ratio_of_P30_A[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Ratio_of_P30_z, Ratio_of_P30_A[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Permeability_x, Permeability_A[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Permeability_y, Permeability_A[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Permeability_z, Permeability_A[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Q_error_x, Q_error_A[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Q_error_y, Q_error_A[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Q_error_z, Q_error_A[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Min_ele_edge_x, min_ele_edge[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Max_ele_edge_x, max_ele_edge[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Mean_ele_edge_x, Mean_ele_edge[0]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Min_ele_edge_y, min_ele_edge[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Max_ele_edge_y, max_ele_edge[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Mean_ele_edge_y, Mean_ele_edge[1]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Min_ele_edge_z, min_ele_edge[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Max_ele_edge_z, max_ele_edge[2]);
    H5file.Append_dataset_to_group(FileKey_mat, groupname, string_Mean_ele_edge_z, Mean_ele_edge[2]);
};

inline void Multi_processes::Matlab_command(string FileKey_m, string FileKey_mat, size_t np, size_t ny, size_t model_no)
{
    std::ofstream oss(FileKey_m, ios::out);
    oss << "clc;\nclose all;\nclear all;\n";
    oss << "s_" << model_no << " = load('" << FileKey_mat << "');\n";

    string L = "L_" + to_string(model_no) + " = s_" + to_string(model_no) + ".Domain_size;\n";
    oss << L;

    string string_P32_total = "P32_total_" + to_string(model_no);
    string string_P30 = "P30_" + to_string(model_no);
    string string_P32_connected_x = "P32_connected_x_" + to_string(model_no);
    string string_P32_connected_y = "P32_connected_y_" + to_string(model_no);
    string string_P32_connected_z = "P32_connected_z_" + to_string(model_no);

    string string_Ratio_of_P32_x = "Ratio_of_P32_x_" + to_string(model_no);
    string string_Ratio_of_P32_y = "Ratio_of_P32_y_" + to_string(model_no);
    string string_Ratio_of_P32_z = "Ratio_of_P32_z_" + to_string(model_no);

    string string_n_I = "n_I_" + to_string(model_no);

    string string_Percolation_probability_x = "Percolation_probability_x_" + to_string(model_no);
    string string_Percolation_probability_y = "Percolation_probability_y_" + to_string(model_no);
    string string_Percolation_probability_z = "Percolation_probability_z_" + to_string(model_no);

    string string_P30_largest_cluster = "P30_largest_cluster_" + to_string(model_no);
    string string_P32_largest_cluster = "P32_largest_cluster_" + to_string(model_no);

    string string_P30_connected_x = "P30_connected_x_" + to_string(model_no);
    string string_P30_connected_y = "P30_connected_y_" + to_string(model_no);
    string string_P30_connected_z = "P30_connected_z_" + to_string(model_no);

    string string_Ratio_of_P30_x = "Ratio_of_P30_x_" + to_string(model_no);
    string string_Ratio_of_P30_y = "Ratio_of_P30_y_" + to_string(model_no);
    string string_Ratio_of_P30_z = "Ratio_of_P30_z_" + to_string(model_no);

    string string_Permeability_x = "Permeability_x_" + to_string(model_no);
    string string_Permeability_y = "Permeability_y_" + to_string(model_no);
    string string_Permeability_z = "Permeability_z_" + to_string(model_no);

    string string_Q_error_x = "Q_error_x_" + to_string(model_no);
    string string_Q_error_y = "Q_error_y_" + to_string(model_no);
    string string_Q_error_z = "Q_error_z_" + to_string(model_no);

    vector<string> String_all = {
        string_P32_total,
        string_P30,
        string_P32_connected_x,
        string_P32_connected_y,
        string_P32_connected_z,
        string_Ratio_of_P32_x,
        string_Ratio_of_P32_y,
        string_Ratio_of_P32_z,
        string_n_I,
        string_Percolation_probability_x,
        string_Percolation_probability_y,
        string_Percolation_probability_z,
        string_P30_largest_cluster,
        string_P32_largest_cluster,
        string_P30_connected_x,
        string_P30_connected_y,
        string_P30_connected_z,
        string_Ratio_of_P30_x,
        string_Ratio_of_P30_y,
        string_Ratio_of_P30_z,
        string_Permeability_x,
        string_Permeability_y,
        string_Permeability_z,
        string_Q_error_x,
        string_Q_error_y,
        string_Q_error_z};

    vector<string> String_matlab_all = {
        "P32_total_",
        "P30_",
        "P32_connected_x_",
        "P32_connected_y_",
        "P32_connected_z_",
        "Ratio_of_P32_x_",
        "Ratio_of_P32_y_",
        "Ratio_of_P32_z_",
        "n_I_",
        "Percolation_probability_x_",
        "Percolation_probability_y_",
        "Percolation_probability_z_",
        "P30_largest_cluster_",
        "P32_largest_cluster_",
        "P30_connected_x_",
        "P30_connected_y_",
        "P30_connected_z_",
        "Ratio_of_P30_x_",
        "Ratio_of_P30_y_",
        "Ratio_of_P30_z_",
        "Permeability_x_",
        "Permeability_y_",
        "Permeability_z_",
        "Q_error_x_",
        "Q_error_y_",
        "Q_error_z_"};

    //cout << String_all.size() << endl;
    //cout << String_matlab_all.size() << endl;

    for (size_t rf = 0; rf < String_all.size() - 3; ++rf)
    {
        if (rf < String_all.size() - 6)
        {
            for (size_t ft = 1; ft <= np; ++ft)
            {
                if (ft == 1)
                    oss << String_all[rf] << " = [mean(s_" << model_no << "." << String_matlab_all[rf] << ft << "); ...\n";
                else if (ft == np)
                    oss << "mean(s_" << model_no << "." << String_matlab_all[rf] << ft << ")];\n";
                else
                    oss << "mean(s_" << model_no << "." << String_matlab_all[rf] << ft << "); ...\n";
            }
            oss << endl;
        }
        else
        {
            for (size_t ft = 1; ft <= np; ++ft)
            {
                oss << "x = find(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << " > 10);\n";
                oss << "s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "(x) = [];\n";

                oss << "if(isempty(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "))\n";
                oss << "\ts_" << model_no << "." << String_matlab_all[rf + 3] << ft << " = [0];\n";
                oss << "end\n\n";

                oss << "s_" << model_no << "." << String_matlab_all[rf] << ft << "(x) = [];\n";
                oss << "if(isempty(s_" << model_no << "." << String_matlab_all[rf] << ft << "))\n";
                oss << "\ts_" << model_no << "." << String_matlab_all[rf] << ft << " = [0];\n";
                oss << "end\n\n";

                //------------------------
                oss << "x = find(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << " == -1);\n";
                oss << "s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "(x) = [];\n";
                oss << "if(isempty(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "))\n";
                oss << "\ts_" << model_no << "." << String_matlab_all[rf + 3] << ft << " = [0];\n";
                oss << "end\n\n";

                oss << "s_" << model_no << "." << String_matlab_all[rf] << ft << "(x) = [];\n";
                oss << "if(isempty(s_" << model_no << "." << String_matlab_all[rf] << ft << "))\n";
                oss << "\ts_" << model_no << "." << String_matlab_all[rf] << ft << " = [0];\n";
                oss << "end\n\n";

                //-------------------isnan
                oss << "x = find(isnan(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "));\n";
                oss << "s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "(x) = [];\n";
                oss << "x = find(isnan(s_" << model_no << "." << String_matlab_all[rf] << ft << "));\n";
                oss << "s_" << model_no << "." << String_matlab_all[rf] << ft << "(x) = [];\n\n";
            }

            oss << endl;

            for (size_t ft = 1; ft <= np; ++ft)
            {
                if (ft == 1)
                    oss << String_all[rf] << " = [mean(s_" << model_no << "." << String_matlab_all[rf] << ft << "); ...\n";
                else if (ft == np)
                    oss << "mean(s_" << model_no << "." << String_matlab_all[rf] << ft << ")];\n";
                else
                    oss << "mean(s_" << model_no << "." << String_matlab_all[rf] << ft << "); ...\n";
            }
            oss << endl;

            for (size_t ft = 1; ft <= np; ++ft)
            {
                if (ft == 1)
                    oss << String_all[rf + 3] << " = [mean(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "); ...\n";
                else if (ft == np)
                    oss << "mean(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << ")];\n";
                else
                    oss << "mean(s_" << model_no << "." << String_matlab_all[rf + 3] << ft << "); ...\n";
            }
            oss << endl;
        }
    }

    oss << "clear s_" << model_no << ";\n";
    oss << "filepath = which('" << FileKey_m << "');\n";
    oss << "filepath = erase(filepath, '" << FileKey_m << "');\n";
    oss << "filepath = [filepath, 'connectivity_modelno_" << model_no << ".mat'];\n";
    oss << "save(filepath);\n";

    oss.close();
};

inline void Multi_processes::Sign_of_finding_pc(string FileKey)
{
    //Writing data
    std::ofstream oss(FileKey, ios::out);
    oss << "P30_c has been found: ";
    oss << Density_c;

    oss.close();
};

inline Vector2d Multi_processes::Re_MHFEM(DFN::Domain dom_222,
                                          const double min_ele_edge,
                                          const double max_ele_edge,
                                          size_t dir_,
                                          double h_in,
                                          double h_out,
                                          size_t Nproc_t_)
{
    string AS = "Fractures_" + To_string_with_width(random_unsigned_integer(1, 1e5), 5) + ".mat";
    dom_222.Matlab_Out_Frac_matfile(AS);

    std::vector<std::vector<Vector3d>> verts__;
    DFN::Load_a_DFN_from_matfile loadmat(AS, verts__);
    DFN::Domain dom_t_;

    Vector6d modelsize__;
    double L_ = this->L;
    modelsize__ << -0.5 * L_, 0.5 * L_, -0.5 * L_, 0.5 * L_, -0.5 * L_, 0.5 * L_;
    dom_t_.Create_whole_model_II(modelsize__, verts__);
    dom_t_.Identify_percolation_clusters();
    dom_t_.Connectivity_analysis();
    dom_t_.Re_identify_intersection_considering_trimmed_frac();
    dom_t_.Identify_percolation_clusters();

    DFN::Mesh_DFN_linear mesh_t_(dom_t_, min_ele_edge * 0.9, max_ele_edge * 0.9, dir_, 1);

    DFN::MHFEM fem_t_{mesh_t_, dom_t_, h_in, h_out, dir_, Nproc_t_};
    Vector2d OSD;
    OSD << fem_t_.Permeability, fem_t_.Q_error;

    const char *CS = AS.c_str();
    std::remove(CS);

    return OSD;
};

}; // namespace DFN
