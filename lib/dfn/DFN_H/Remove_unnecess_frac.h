#pragma once
#include "Domain_WL.h"

namespace DFN
{
class Remove_unnecess_frac
{
public:
    std::vector<std::vector<int>> Listofclusters_remove_fracs;

public:
    Remove_unnecess_frac(DFN::Domain dom);
};

inline Remove_unnecess_frac::Remove_unnecess_frac(DFN::Domain dom)
{
    Listofclusters_remove_fracs.resize(dom.Percolation_cluster[2].size());

    for (size_t i = 0; i < dom.Percolation_cluster[2].size(); ++i)
    {
        Listofclusters_remove_fracs[i].resize(dom.Listofclusters[dom.Percolation_cluster[2][i]].size());
        for (size_t j = 0; j < dom.Listofclusters[dom.Percolation_cluster[2][i]].size(); ++j)
        {
            Listofclusters_remove_fracs[i][j] = dom.Listofclusters[dom.Percolation_cluster[2][i]][j];
        }
    }
/*
    for (size_t i = 0; i < Listofclusters_remove_fracs.size(); ++i)
    {
        for (size_t j = 0; j < Listofclusters_remove_fracs[i].size(); ++j)
        {
            size_t fracID = Listofclusters_remove_fracs[i][j];
            size_t Intersect_num = 0;

            if (dom.Fractures[fracID].If_intersect_surfaces[0] == 1)
                Intersect_num++;

            if (dom.Fractures[fracID].If_intersect_surfaces[1] == 1)
                Intersect_num++;

            Intersect_num += dom.Fractures[fracID].Intersect_other_frac_after_trim.size();

            if (Intersect_num < 2)
            {
                Listofclusters_remove_fracs[i][j] = -1;
            }
        }
    }
*/
};
}; // namespace DFN