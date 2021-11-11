#pragma once
#include "mat.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace DFN
{
class MATLAB_DATA_API
{
public:
    MATLAB_DATA_API();
    void Write_mat(string FileKey_mat, string mode, size_t NUM_eles, size_t rows, size_t cols, vector<double> data_, /*column major data*/
              string field_name);                                                                               // only for at most 2D matrix
};

inline MATLAB_DATA_API::MATLAB_DATA_API(){

};

inline void MATLAB_DATA_API::Write_mat(string FileKey_mat, string mode, size_t NUM_eles, size_t rows, size_t cols, vector<double> data_, /*column major data*/
                                  string field_name)
{
    const char *filename = FileKey_mat.c_str();
    MATFile *pMatFile;
    const char *mode_ = mode.c_str();
    pMatFile = matOpen(filename, mode_);

    if (!pMatFile)
        throw Error_throw_pause("cannot create mat file\n");

    double pData1[NUM_eles] = {};

    mxArray *pMxArray1;
    pMxArray1 = mxCreateDoubleMatrix(rows, cols, mxREAL);

    if (!pMxArray1 || !pData1)
    {
        matClose(pMatFile);
        throw Error_throw_pause("cannot create pMxArray or pData\n");
    }

    for (size_t j = 0; j < NUM_eles; ++j)
        pData1[j] = data_[j];

    memcpy((void *)(mxGetPr(pMxArray1)), (void *)pData1, sizeof(pData1));

    const char *field_ = field_name.c_str();

    matPutVariable(pMatFile, field_, pMxArray1);

    mxDestroyArray(pMxArray1);

    matClose(pMatFile);
};

} // namespace DFN