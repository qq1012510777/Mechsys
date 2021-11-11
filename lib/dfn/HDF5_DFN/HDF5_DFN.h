#pragma once

#include "../Error_throw/Error_throw.h"
#include "H5Cpp.h"
#include "hdf5.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace H5;

namespace DFN
{
class HDF5_DFN
{
public:
    HDF5_DFN();
    HDF5_DFN(string filename);
    void Write_H5(string filename, string groupname, string datasetname, vector<double> data);
    void Write_H5(string filename, string datasetname, string data);
    void Write_H5(string filename, string datasetname, double data);
    void Write_H5(string filename, string groupname, vector<string> datasetname, vector<vector<double>> data);
    void Overwrite(string filename, string datasetname, double data);
    vector<double> Read_H5(string filename, string groupname, string datasetname);
    vector<double> Read_H5(string filename, string datasetname);
    string Read_H5_text(string filename, string datasetname);
    void Append_dataset_to_group(string filename, string groupname, string datasetname, vector<double> data);
};

HDF5_DFN::HDF5_DFN(){

};

HDF5_DFN::HDF5_DFN(string filename)
{
    H5File file(filename, H5F_ACC_TRUNC);
    file.close();
};

///1 double
void HDF5_DFN::Write_H5(string filename, string groupname, string datasetname, vector<double> data)
{
    size_t dim = data.size();

    H5File file(filename, H5F_ACC_RDWR);
    Group group(file.createGroup(groupname));
    hsize_t dims[1]; // dataset dimensions for each rank
    dims[0] = dim;
    // Create the dataspace for a dataset first.
    DataSpace dataspace(1, dims);
    //IntType datatype(PredType::NATIVE_FLOAT);
    //datatype.setOrder(H5T_ORDER_LE);
    // Create the dataset under group with specified dataspace.
    DataSet dataset = group.createDataSet(datasetname, PredType::NATIVE_DOUBLE, dataspace);
    double *buffer = new double[dim]();
    if (buffer == NULL)
    {
        string AS = "in Write_H5(), error buff\n";
        throw Error_throw_pause(AS);
    }
    for (size_t i = 0; i < dim; ++i)
        buffer[i] = data[i];
    dataset.write(buffer, PredType::NATIVE_DOUBLE);

    delete[] buffer;
    buffer = NULL;

    group.close();
    file.close();
};

///2 string
void HDF5_DFN::Write_H5(string filename, string datasetname, string data)
{
    H5File file(filename, H5F_ACC_RDWR);

    H5::StrType h5stringType(H5::PredType::C_S1, data.length() + 1); // + 1 for trailing zero
    H5::DataSet ds = file.createDataSet(datasetname, h5stringType, H5::DataSpace(H5S_SCALAR));
    ds.write(data, h5stringType);

    file.close();
};

void HDF5_DFN::Overwrite(string filename, string datasetname, double data)
{
    H5File file(filename, H5F_ACC_RDWR); //The hdf5 c++ object.
    string channelName = "/" + datasetname;
    int result = H5Ldelete(file.getId(), channelName.data(), H5P_DEFAULT);
    result++;
    file.close();

    this->Write_H5(filename, datasetname, data);
}

///3 one double value
void HDF5_DFN::Write_H5(string filename, string datasetname, double data)
{
    //size_t RANK = 1;

    H5File file(filename, H5F_ACC_RDWR);

    hsize_t dims[1]; // dataset dimensions for each rank
    dims[0] = 1;

    DataSpace dataspace(1, dims);

    DataSet dataset = file.createDataSet(datasetname, PredType::NATIVE_DOUBLE, dataspace);
    double buffer[1] = {data};
    dataset.write(buffer, PredType::NATIVE_DOUBLE);

    file.close();
};

void HDF5_DFN::Write_H5(string filename, string groupname, vector<string> datasetname, vector<vector<double>> data)
{
    H5File file(filename, H5F_ACC_RDWR);
    Group group(file.createGroup(groupname));

    for (size_t k = 0; k < datasetname.size(); ++k)
    {
        hsize_t dims[1]; // dataset dimensions for each rank
        dims[0] = data[k].size();
        // Create the dataspace for a dataset first.
        DataSpace dataspace(1, dims);
        //IntType datatype(PredType::NATIVE_FLOAT);
        //datatype.setOrder(H5T_ORDER_LE);
        // Create the dataset under group with specified dataspace.

        DataSet dataset = group.createDataSet(datasetname[k], PredType::NATIVE_DOUBLE, dataspace);

        double *buffer = new double[data[k].size()]();
        if (buffer == NULL)
        {
            string AS = "in Write_H5(), error buff\n";
            throw Error_throw_pause(AS);
        }

        for (size_t i = 0; i < data[k].size(); ++i)
            buffer[i] = data[k][i];
        dataset.write(buffer, PredType::NATIVE_DOUBLE);

        delete[] buffer;
        buffer = NULL;
    }

    group.close();
    file.close();
};

inline vector<double> HDF5_DFN::Read_H5(string filename, string groupname, string datasetname)
{
    H5File file(filename, H5F_ACC_RDWR);

    Group group = file.openGroup(groupname);

    DataSet dataset = group.openDataSet(datasetname);

    DataSpace filespace = dataset.getSpace();

    int rank = filespace.getSimpleExtentNdims();

    if (rank != 1)
    {
        string AS = "Sorry, this function is only for 1D data!\n";
        throw Error_throw_pause(AS);
    }

    hsize_t dims[rank];

    rank = filespace.getSimpleExtentDims(dims);

    DataSpace myspace(rank, dims);

    double *buffer = new double[dims[0]]();

    dataset.read(buffer, PredType::NATIVE_DOUBLE, myspace, filespace);

    vector<double> AY(dims[0]);

    for (size_t i = 0; i < dims[0]; ++i)
        AY[i] = buffer[i];

    delete[] buffer;
    buffer = NULL;
    group.close();
    file.close();

    return AY;
};

inline vector<double> HDF5_DFN::Read_H5(string filename, string datasetname)
{
    H5File file(filename, H5F_ACC_RDWR);

    DataSet dataset = file.openDataSet(datasetname);

    DataSpace filespace = dataset.getSpace();

    int rank = filespace.getSimpleExtentNdims();

    if (rank != 1)
    {
        string AS = "Sorry, this function is only for 1D data!\n";
        throw Error_throw_pause(AS);
    }

    hsize_t dims[rank];

    rank = filespace.getSimpleExtentDims(dims);

    DataSpace myspace(rank, dims);

    double *buffer = new double[dims[0]]();

    dataset.read(buffer, PredType::NATIVE_DOUBLE, myspace, filespace);

    vector<double> AY(dims[0]);

    for (size_t i = 0; i < dims[0]; ++i)
        AY[i] = buffer[i];

    delete[] buffer;
    buffer = NULL;
    file.close();

    return AY;
};

inline string HDF5_DFN::Read_H5_text(string filename, string datasetname)
{
    H5File file(filename, H5F_ACC_RDWR);

    // open dataset, get data-type
    DataSet dataset = file.openDataSet(datasetname);
    DataSpace dataspace = dataset.getSpace();
    StrType datatype = dataset.getStrType();

    // allocate output
    string data;

    // read output
    dataset.read(data, datatype, dataspace);

    return data;
};

inline void HDF5_DFN::Append_dataset_to_group(string filename, string groupname, string datasetname, vector<double> data)
{
    H5File file(filename, H5F_ACC_RDWR);

    Group group = file.openGroup(groupname);

    hsize_t dims[1]; // dataset dimensions for each rank
    dims[0] = data.size();

    DataSpace dataspace(1, dims);

    DataSet dataset = group.createDataSet(datasetname, PredType::NATIVE_DOUBLE, dataspace);

    double *buffer = new double[data.size()]();

    if (buffer == NULL)
    {
        string AS = "in Append_dataset_to_group(), error buff\n";
        throw Error_throw_pause(AS);
    }

    for (size_t i = 0; i < data.size(); ++i)
        buffer[i] = data[i];

    dataset.write(buffer, PredType::NATIVE_DOUBLE);

    delete[] buffer;
    buffer = NULL;

    group.close();
    file.close();
};

}; // namespace DFN