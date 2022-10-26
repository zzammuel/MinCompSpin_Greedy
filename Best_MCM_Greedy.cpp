// To compile: g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp info_quant.cpp best_basis.cpp Basis_Choice.cpp P_s.cpp metropolis.cpp
// To run: time ./a.out
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <sstream>
#include <bitset>
#include <map>

using namespace std;

/******************************************************************************/
/*******************************    CONSTANTS     *****************************/
/******************************************************************************/
#include "data.h"

/******************************************************************************/
/*****************************    PRINT MCM   *********************************/
/******************************************************************************/
//unsigned int A = 5;   // A = number of parts
//uniform_int_distribution<int> uni(0., A - 1);

void Print_Partition(map<unsigned int, __int128_t> partition)
{
    map<unsigned int, __int128_t>::iterator it;

    for (it = partition.begin(); it != partition.end(); it++)
    {
        bitset<n> hi{ static_cast<unsigned long long>((*it).second >> 64) },
            lo{ static_cast<unsigned long long>((*it).second) },
            bits{ (hi << 64) | lo };
        cout << (*it).first << "\t " << bits << endl;
    }
    cout << endl;
}

/******************************************************************************/
/***************    FIND THE BEST MCM: GREEDY PROCEDURE    ********************/
/******************************************************************************/

double LogE_SubCM(map<__int128_t, unsigned int > Kset, __int128_t Ai, unsigned int N, bool print_bool = false);


map<unsigned int, __int128_t> Matrix(map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    // Create initial independent model:
    map<unsigned int, __int128_t> Best_MCM_Partition;  __int128_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_Part;

    cout << "**** Filling in the nxn-triangular Matrix: ****" << endl;
    double* logE_Mat = (double*)malloc(r * (r - 1) / 2 * sizeof(double));

    for (int i = 0; i < r - 1; i++)
    {
        for (int j = i + 1; j < r; j++)
        {
            logE_Mat[j * (j - 1) / 2 + i] = 0;
        }
    }

    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = Op;
        Best_MCM_LogE_Part[i] = LogE_SubCM(Kset, Op, N);
        Op = Op << 1;
    }

    // buffers:
    map<unsigned int, __int128_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, __int128_t>::iterator it1, it2, it1_end, it2_start, it;
    int counter = 0;
    bool stop = false;

    int i_keep, i_erase, i_mat; __int128_t new_part;   // define variables for merging

    int k = 0, iteration = 0;
    
    cout << "**** Start Merging: ****" << endl;
    // ********* HIERARCHICAL MERGING ****************
    while (Best_MCM_Partition.size() > 1 && (!stop))
    {
        counter = 0;
        Diff_LogE_best = 0;
        LogE_merged = 0, LogE_unmerged = 0;

        // ****** Find the best pair to merge: *********
        //cout << endl << "******** Test combinations: " << endl;
        it1_end = Best_MCM_Partition.end();   it1_end--;
        for (it1 = Best_MCM_Partition.begin(); it1 != it1_end; it1++)
        {
            it2_start = it1;    it2_start++;
            for (it2 = it2_start; it2 != Best_MCM_Partition.end(); it2++)
            {
                i_mat = (*it2).first * ((*it2).first - 1) / 2 + (*it1).first;

                //cout << counter << "\t" << bitset<n>((*it1).second) << "\t " << bitset<n>((*it2).second);

                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

                if (logE_Mat[i_mat] == 0)
                {
                    logE_Mat[i_mat] = LogE_SubCM(Kset, (*it1).second + (*it2).second, N);
                }
                //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged);

                if ((logE_Mat[i_mat] - LogE_unmerged) > Diff_LogE_best)
                {
                    Buffer_Community = MCM_unmerged;
                    Diff_LogE_best = logE_Mat[i_mat] - LogE_unmerged;
                    LogE_merged_best = logE_Mat[i_mat];
                }
                //cout << "\t " << Diff_LogE_best << endl;
                //counter++;
                MCM_unmerged.clear();
            }
        }

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; continue; }

        // ********* PERFORM MERGING *****:
        it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
        it++;  i_erase = (*it).first; //  new_part += (*it).second;

        for (it1 = Best_MCM_Partition.begin(); it1 != Best_MCM_Partition.end(); it1++)
        {
            if ((*it1).first == i_keep) { continue; }
            if ((*it1).first > i_keep)
            {
                logE_Mat[(*it1).first * ((*it1).first - 1) / 2 + i_keep] = 0;
            }
            else
            {
                logE_Mat[i_keep * (i_keep - 1) / 2 + (*it1).first] = 0;
            }

        }

        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        iteration++;
        cout << "Just done with iteration " << iteration << endl;
    } // End While

    return Best_MCM_Partition;
}


map<unsigned int, __int128_t> Matrix2(map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    // Create initial independent model:
    map<unsigned int, __int128_t> Best_MCM_Partition;  __int128_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_Part;

    double* logE_Mat = (double*)malloc(r * (r - 1) / 2 * sizeof(double));

    for (int i = 0; i < r - 1; i++)
    {
        for (int j = i + 1; j < r; j++)
        {
            logE_Mat[j * (j - 1) / 2 + i] = 0;
        }
    }

    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = Op;
        Best_MCM_LogE_Part[i] = LogE_SubCM(Kset, Op, N);
        Op = Op << 1;
    }

    // buffers:
    map<unsigned int, __int128_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, __int128_t>::iterator it1, it2, it1_end, it2_start, it;
    int counter = 0;
    bool stop = false;

    int i_keep, i_erase, i_mat, i_mat2; __int128_t new_part;   // define variables for merging

    int k = 0, iteration = 0;

    // ********* HIERARCHICAL MERGING ****************
    while (Best_MCM_Partition.size() > 1 && (!stop))
    {
        Diff_LogE_best = 0;
        LogE_merged = 0, LogE_unmerged = 0;

        // ****** Find the best pair to merge: *********
        it1_end = Best_MCM_Partition.end();   it1_end--;
        for (it1 = Best_MCM_Partition.begin(); it1 != it1_end; it1++)
        {
            it2_start = it1;    it2_start++;
            for (it2 = it2_start; it2 != Best_MCM_Partition.end(); it2++)
            {
                i_mat = (*it2).first * ((*it2).first - 1) / 2 + (*it1).first;

                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

                if (logE_Mat[i_mat] == 0)
                {
                    logE_Mat[i_mat] = LogE_SubCM(Kset, (*it1).second + (*it2).second, N) - LogE_unmerged;
                }

                if ((logE_Mat[i_mat]) > Diff_LogE_best)
                {
                    Buffer_Community = MCM_unmerged;
                    Diff_LogE_best = logE_Mat[i_mat];
                    LogE_merged_best = logE_Mat[i_mat] + LogE_unmerged;
                }
                MCM_unmerged.clear();
            }
        }

        // ********* PERFORM MERGING *****:
        it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
        it++;  i_erase = (*it).first; //  new_part += (*it).second;

        for (it1 = Best_MCM_Partition.begin(); it1 != Best_MCM_Partition.end(); it1++)
        {
            if (((*it1).first == i_keep) || ((*it1).first == i_erase))
            {
                continue;
            }

            if ((*it1).first > i_keep)
            {
                i_mat = (*it1).first * ((*it1).first - 1) / 2 + i_keep;
            }
            else
            {
                i_mat = i_keep * (i_keep - 1) / 2 + (*it1).first;
            }

            if ((*it1).first > i_erase)
            {
                i_mat2 = (*it1).first * ((*it1).first - 1) / 2 + i_erase;
            }
            else
            {
                i_mat2 = i_erase * (i_erase - 1) / 2 + (*it1).first;
            }

            if ((logE_Mat[i_mat] >= 0) || (logE_Mat[i_mat2] >= 0))
            {
                logE_Mat[i_mat] = 0;
            }
            else { counter++; }
        }

        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; }
        iteration++;
    } // End While

    return Best_MCM_Partition;
}
