#pragma once
#include "../Error_throw/Error_throw.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <map>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{
class Triplet_Eigen
{
public:
    vector<Triplet<double>> tripletList;
    size_t size_triplet;
    size_t push_time_counter = 0;
    map<pair<size_t, size_t>, size_t> glob_idx;
    size_t global_IDX_counter = 0;

public:
    Triplet_Eigen();
    Triplet_Eigen(const size_t size_tp);

    void Fill_triplet(const size_t glob_id,
                      const size_t i,
                      const size_t j,
                      const double value);

    void Fill_triplet_duplicate(const size_t i,
                                const size_t j,
                                const double value);
    void Clear_triplet();
};

inline Triplet_Eigen::Triplet_Eigen(){};

inline Triplet_Eigen::Triplet_Eigen(const size_t size_tp)
{
    this->tripletList.resize(size_tp);
    this->size_triplet = size_tp;
};

inline void Triplet_Eigen::Fill_triplet(const size_t glob_id,
                                        const size_t i,
                                        const size_t j,
                                        const double value)
{
    if (glob_id < this->size_triplet)
        tripletList[glob_id] = Triplet<double>(j, i, value);
    else
    {
        string AS = "the index exceeds the range of the vector 'Fill_triplet', in class of 'Triplet_Eigen'!\n";
        throw Error_throw_pause(AS);
    }
};

inline void Triplet_Eigen::Fill_triplet_duplicate(const size_t i,
                                                  const size_t j,
                                                  const double value)
{
    size_t glob_id = 0;

    map<pair<size_t, size_t>, size_t>::iterator its;
    its = this->glob_idx.find(std::make_pair(i, j));

    if (its != this->glob_idx.end()) // duplicate idx
        glob_id = its->second;
    else // not a duplicate point
    {
        glob_idx.insert(std::make_pair(make_pair(i, j), this->global_IDX_counter));
        glob_id = this->global_IDX_counter;
        this->global_IDX_counter++;
    };

    if (glob_id < this->size_triplet)
    {
        double value_s = tripletList[glob_id].value() + value;
        tripletList[glob_id] = Triplet<double>(j, i, value_s);
    }
    else
    {
        string AS = "the index exceeds the range of the vector 'Fill_triplet_duplicate', in class of 'Triplet_Eigen'!\n";
        throw Error_throw_pause(AS);
    }
};

inline void Triplet_Eigen::Clear_triplet()
{
    this->tripletList.clear();
};

}; // namespace DFN
