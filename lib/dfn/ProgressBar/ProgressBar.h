#pragma once

#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

namespace DFN
{
class ProgressBar
{
public:
    size_t global_counter; // for parallelization
    vector<double> Spacing_counter;
    size_t Lo_counter = 0;

public:
    ProgressBar();
    void Rep_prog_serially_for_terminal(size_t i /*step*/, size_t Total_num, string As);
    void Rep_prog_serially_for_supercomputer(size_t i /*step*/, size_t spacing, size_t Total_num, string As);
    void Rep_prog_for_paral(size_t Total_num, size_t spacing, string As);
    ~ProgressBar();
};

//--------------------
ProgressBar::ProgressBar()
{
    global_counter = 0;
    Spacing_counter.resize(1);
    Spacing_counter[0] = 0;
};

inline void ProgressBar::Rep_prog_serially_for_terminal(size_t i /*step*/, size_t Total_num, string As)
{
    cout << "\r" << As;

    printf(" [%6.2lf%%]:", i * 100.0 / (Total_num - 1));

    int show_num = i * 20 / Total_num;

    for (int j = 1; j <= show_num; j++)
        cout << "█";
};

inline void ProgressBar::Rep_prog_serially_for_supercomputer(size_t i /*step*/, size_t spacing, size_t Total_num, string As)
{
    double prog = i * 1.0 / (Total_num * 1.0);

    if (prog >= Spacing_counter[Lo_counter])
    {
        Lo_counter++;
        Spacing_counter.push_back(Spacing_counter[Lo_counter - 1] + (1.0 / spacing));
        cout << As;
        printf(" [%6.2lf%%]:", prog * 100.0);

        int show_num = i * 20.0 / Total_num;

        for (int j = 1; j <= show_num; j++)
            cout << "█";
        //cout << Spacing_counter[Lo_counter - 1];
        cout << endl;
    }
};

inline void ProgressBar::Rep_prog_for_paral(size_t Total_num, size_t spacing, string As)
{
#pragma omp critical
    {
        global_counter++;

        double prog = global_counter * 1.0 / (Total_num * 1.0);

        //cout << "prog " << prog << ", " << Spacing_counter[Lo_counter] << endl;

        if (prog >= Spacing_counter[Lo_counter])
        {
            Lo_counter++;
            Spacing_counter.push_back(Spacing_counter[Lo_counter - 1] + (1.0 / spacing));

            cout << As << ": ";

            size_t i = global_counter;

            printf("[%6.2f%%]:", i * 100.0 / (Total_num));

            int show_num = i * 20 / Total_num;

            for (int j = 1; j <= show_num; j++)
                cout << "█";
            cout << endl;
        }
    };
};

ProgressBar::~ProgressBar()
{
}

} // namespace DFN