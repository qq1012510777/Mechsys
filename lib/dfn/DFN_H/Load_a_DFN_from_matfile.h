#pragma once
#include "../Error_throw/Error_throw.h"
#include "Eigen/Dense"
#include "mat.h"
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{

class Load_a_DFN_from_matfile
{
public:
    Load_a_DFN_from_matfile(string FileKey_mat, std::vector<std::vector<Vector3d>> &verts);
};

Load_a_DFN_from_matfile::Load_a_DFN_from_matfile(string FileKey_mat, std::vector<std::vector<Vector3d>> &verts)
{
    const char *filename = FileKey_mat.c_str();
    MATFile *pMatFile;
    pMatFile = matOpen(filename, "r");

    if (!pMatFile)
    {
        throw Error_throw_ignore("cannot create mat file in class Domain\n");
    }

    mxArray *pMxArray = NULL;
    pMxArray = matGetVariable(pMatFile, "Num_fracs");
    double *getdata;
    size_t Num_fractures = 0;

    if (!pMxArray)
    {
        cout << "Cannot find fractures matfile!\n";
        throw Error_throw_pause("Cannot find fractures matfile!\n");
    }
    else
    {
        getdata = (double *)mxGetData(pMxArray);
        Num_fractures = getdata[0];
    }

    mxFree(getdata);
    //cout << "Input Num_fractures: " << Num_fractures << endl;
    verts.resize(Num_fractures);
    if (Num_fractures != 0)
        for (size_t i = 0; i < Num_fractures; ++i)
        {
            string ft = to_string(i + 1);
            string Fracx = "Frac_" + ft + "_x";
            const char *Fracx_s = Fracx.c_str();

            string Fracy = "Frac_" + ft + "_y";
            const char *Fracy_s = Fracy.c_str();

            string Fracz = "Frac_" + ft + "_z";
            const char *Fracz_s = Fracz.c_str();

            mxArray *pMxArray1 = NULL;
            mxArray *pMxArray2 = NULL;
            mxArray *pMxArray3 = NULL;
            pMxArray1 = matGetVariable(pMatFile, Fracx_s);
            pMxArray2 = matGetVariable(pMatFile, Fracy_s);
            pMxArray3 = matGetVariable(pMatFile, Fracz_s);
            double *getdata1, *getdata2, *getdata3;

            getdata1 = (double *)mxGetData(pMxArray1);
            getdata2 = (double *)mxGetData(pMxArray2);
            getdata3 = (double *)mxGetData(pMxArray3);

            std::vector<Vector3d> verts_kkk(4);
            for (size_t j = 0; j < 4; ++j)
            {
                verts_kkk[j] << getdata1[j], getdata2[j], getdata3[j];
            }
            verts[i] = verts_kkk;
            mxFree(getdata1);
            mxFree(getdata2);
            mxFree(getdata3);
        }

    matClose(pMatFile);
};

}; // namespace DFN