#pragma once
#include "../DFN_H/Domain_WL.h"
#include "../Error_throw/Error_throw.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

namespace DFN
{
class Remove_fractures
{
public:
    Remove_fractures();

public:
    bool Remove_dead_fractures(std::vector<size_t> &One_cluster, DFN::Domain dom, size_t dir);
};

inline Remove_fractures::Remove_fractures(){

};

inline bool Remove_fractures::Remove_dead_fractures(std::vector<size_t> &One_cluster, DFN::Domain dom, size_t dir)
{
    size_t inlet_id = 0, outlet_id = 0;
    if (dir == 0)
    {
        inlet_id = 4;
        outlet_id = 5;
    }
    else if (dir == 1)
    {
        inlet_id = 2;
        outlet_id = 3;
    }
    else if (dir == 2)
    {
        inlet_id = 0;
        outlet_id = 1;
    };

    size_t frac_max = *std::max_element(One_cluster.begin(), One_cluster.end());

    Eigen::SparseMatrix<size_t> adjacent_matrix(frac_max + 3, frac_max + 3);
    adjacent_matrix.reserve(VectorXi::Constant(frac_max + 3, 5));

    for (size_t i = 0; i < One_cluster.size(); ++i)
    {
        size_t FracTag = One_cluster[i];

        if (dom.Fractures[FracTag].If_intersect_surfaces[inlet_id] == 1)
        {
            adjacent_matrix.insert(FracTag, frac_max + 1) = 1;
            adjacent_matrix.insert(frac_max + 1, FracTag) = 1; // the last two
        }

        if (dom.Fractures[FracTag].If_intersect_surfaces[outlet_id] == 1)
        {
            adjacent_matrix.insert(FracTag, frac_max + 2) = 1;
            adjacent_matrix.insert(frac_max + 2, FracTag) = 1; // the last one
        }
    }

    for (size_t i = 0; i < One_cluster.size() - 1; ++i)
    {
        size_t FracTag_1 = One_cluster[i];

        for (size_t j = i + 1; j < One_cluster.size(); ++j)
        {
            size_t FracTag_2 = One_cluster[j];

            //-------
            size_t FracTag_A = FracTag_1 < FracTag_2 ? FracTag_1 : FracTag_2;
            size_t FracTag_B = FracTag_1 > FracTag_2 ? FracTag_1 : FracTag_2;

            std::map<std::pair<size_t, size_t>, std::pair<Vector3d, Vector3d>>::iterator ity;

            ity = dom.Intersections.find(std::make_pair(FracTag_A, FracTag_B));

            if (ity != dom.Intersections.end())
            {
                adjacent_matrix.insert(FracTag_A, FracTag_B) = 1;
                adjacent_matrix.insert(FracTag_B, FracTag_A) = 1;
            }
        }
    }

    adjacent_matrix.makeCompressed();

    bool tees = false;

    for (size_t i = 1; i <= frac_max; ++i)
    {
        SparseVector<size_t> sub_mat = adjacent_matrix.innerVector(i);

        if (sub_mat.nonZeros() == 1) // only connect to one frac
        {
            std::vector<size_t>::iterator itr = std::find(One_cluster.begin(), One_cluster.end(), i);
            size_t local_index = std::distance(One_cluster.begin(), itr);

            tees = true;

            One_cluster.erase(One_cluster.begin() + local_index);
            //cout << i << endl;
            break;
        }
    }

    return tees;
};

}; // namespace DFN