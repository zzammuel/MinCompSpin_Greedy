//g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Basis_Choice.cpp RandomPart.cpp 
// Best_MCM.cpp 
//time ./a.out
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <numeric>
#include <bitset>
#include <map>
#include <cmath>       /* tgamma */
#include <random>
#include <algorithm>
#include <queue> /* priority queue (max heap) */
#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"
#include "library.h"

using namespace std::chrono;
/******************************************************************************/
/**************************     RANDOM PARTITIONS    **************************/
/******************************************************************************/
//std::mt19937 gen;       //only for mt19937    
std::default_random_engine gen;
void initialise_generator()
{
    int seed = (unsigned)time(NULL);
    //srand48(seed);      //for drand48
    gen.seed(seed);     //for mt19937
}

//map<unsigned int, uint64_t> Random_Partition(unsigned int A, mt19937 gen);

unsigned int A = 5;   // A = number of parts
uniform_int_distribution<int> uni(0., A - 1);
uniform_real_distribution<double> uni_r(0.0, 1.0);

map<unsigned int, uint64_t> Random_Partition()
{
    //  uniform_int_distribution<int> uni(0.,A-1);
    map<unsigned int, uint64_t> Partition_rand;
    uint64_t si = 0;
    unsigned int part = 0;

    for (unsigned int a = 0; a < A; a++)
    {
        Partition_rand[a] = si;
    }

    si = 1;

    for (int i = 0; i < n; i++)
    {
        part = uni(gen);
        Partition_rand[part] += si;
        si = (si << 1);
    }
    /*
      for(map<unsigned int, uint64_t>::iterator it = Partition_rand.begin(); it != Partition_rand.end(); it++)
      {
        cout << (*it).second << "\t " << bitset<n>((*it).second) << endl;
      }

      cout << "Is this partition correct?  " << ((check_partition(Partition_rand))?"Yes":"No") << endl;
    */

    return Partition_rand;
}


/********************************************************************/
/**************************    PRINT INFO    ************************/
/******    ON SUCCESSIVE INDEPENDENT MODELS IN THE NEW BASIS   ******/
/********************************************************************/
/*
void Full_Indep_Model(map<uint32_t, unsigned int> Kset, unsigned int N)
{
  map<uint32_t, uint32_t> Partition_Indep;  uint32_t Op = 1;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_Indep[i] = Op;
    cout << "Add Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;
    Op = Op << 1;
  }
  Partition_Indep.clear();
}*/

// *** map<uint32_t, uint32_t>   --> .first = i = index    --> .second = a[i] = number of element in the part
map<uint32_t, uint32_t> Convert_Partition_forMCM(uint32_t* a, unsigned int r = n)
{
    map<uint32_t, uint32_t> Partition;
    uint32_t element = 1;

    for (int i = r - 1; i >= 0; i--)  // read element from last to first
    {
        Partition[(a[i])] += element;
        element = element << 1;      //cout << a[i] << "\t :";   Print_Partition_Converted(Partition);
    }

    //  cout << "Convert, " << Partition.size() << " parts: \t " ;
    //  bool ok = check_partition(Partition);
    //  cout << " \t --> Partition Ok? " << ok << endl << endl;

    return Partition;
}

double return_max(map<unsigned int, map<unsigned int, double>> A, int& i, int& j)
{
    double maximum = 0;
    for (auto const& ent1 : A)
    {
        for (auto const& ent2 : ent1.second)
        {
            if (ent2.second > maximum)
            {
                maximum = ent2.second;
                i = ent1.first;
                j = ent2.first;
            }
        }
    }
    return maximum;
}

void print_keys(map<unsigned int, map<unsigned int, double>> A)
{
    for (auto const& ent1 : A)
    {
        for (auto const& ent2 : ent1.second)
        {
            cout << "(" << ent1.first << ", " << ent2.first << "), ";
        }
        cout << endl;
    }
}

void print_values(map<unsigned int, map<unsigned int, double>> A)
{
    for (auto const& ent1 : A)
    {
        for (auto const& ent2 : ent1.second)
        {
            cout << ent2.second << ", ";
        }
        cout << endl;
    }
}

//double LogE_MCM(map<uint64_t, unsigned int> Kset, map<unsigned int, uint64_t> Partition, unsigned int N, bool print_bool = false);

double LogE_MCM_converted(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int* a, unsigned int r = n)
{
    map<unsigned int, uint64_t> Partition;
    uint64_t element = 1;

    for (int i = r - 1; i >= 0; i--)  // read element from last to first
    {
        Partition[(a[i])] += element;
        element = element << 1;      //cout << a[i] << "\t :";   Print_Partition_Converted(Partition);
    }

    double LogE = LogE_MCM(Kset, Partition, N);
    Partition.clear();

    return LogE;
}

void Print_Partition(uint32_t* a, unsigned int r = n)
{
    //for (int i=0; i<n-r; i++)
    //{   cout << "x" << " ";   }
    for (int i = 0; i < r; i++)
    {
        cout << a[i] << " ";
    }
    cout << endl;
}

// Original algorithm
void Original(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    // Create initial independent model:
    map<unsigned int, uint64_t> Best_MCM_Partition;  uint64_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_Part;

    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = Op;
        Best_MCM_LogE_Part[i] = LogE_SubCM(Kset, Op, N);
        Op = Op << 1;
    }

    // buffers:
    map<unsigned int, uint64_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, uint64_t>::iterator it1, it2, it1_end, it2_start, it;
    int counter = 0;
    bool stop = false;

    int i_keep, i_erase; uint64_t new_part;   // define variables for merging

    int k = 0, iteration = 0;

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
                //cout << counter << "\t" << bitset<n>((*it1).second) << "\t " << bitset<n>((*it2).second);

                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

                LogE_merged = LogE_SubCM(Kset, (*it1).second + (*it2).second, N);

                //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged);

                if ((LogE_merged - LogE_unmerged) > Diff_LogE_best)
                {
                    Buffer_Community = MCM_unmerged;
                    Diff_LogE_best = LogE_merged - LogE_unmerged;
                    LogE_merged_best = LogE_merged;
                }
                //cout << "\t " << Diff_LogE_best << endl;
                //counter++;
                MCM_unmerged.clear();
            }
        }
        // ********* PERFORM MERGING *****:
        it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
        it++;  i_erase = (*it).first; //  new_part += (*it).second;

        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; }
        iteration++;
    } // End While

    // ********* FINAL PARTITION **************:
    cout << endl << "******** FINAL PARTITION -- iteration " << iteration << ":" << endl << endl;

    LogE_tot = 0;
    for (it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
    {
        cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
        LogE_tot += Best_MCM_LogE_Part[(*it).first];
    }

    cout << endl << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
    cout << "\t LogE = " << LogE_tot << endl;   // LogE_MCM(Kset, Best_MCM_Partition, N) 
}

// Regular merging with no double log-evidence computations.
map<unsigned int, uint64_t> Matrix(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    // Create initial independent model:
    map<unsigned int, uint64_t> Best_MCM_Partition;  uint64_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_Part;

    double** Diff_Mat = (double**)malloc(r * sizeof(double*));

    for (int i = 0; i < r; i++)
    {
        Diff_Mat[i] = (double*)malloc(r * sizeof(double));
        for (int j = 0; j < r; j++)
        {
            Diff_Mat[i][j] = 0;
        }
    }

    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = Op;
        Best_MCM_LogE_Part[i] = LogE_SubCM(Kset, Op, N);
        Op = Op << 1;
    }

    // buffers:
    map<unsigned int, uint64_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, uint64_t>::iterator it1, it2, it1_end, it2_start, it;
    int counter = 0;
    bool stop = false;

    int i_keep, i_erase; uint64_t new_part;   // define variables for merging

    int k = 0, iteration = 0;

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

                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

                if (Diff_Mat[(*it1).first][(*it2).first] == 0)
                {
                    Diff_Mat[(*it1).first][(*it2).first] = LogE_SubCM(Kset, (*it1).second + (*it2).second, N);
                }

                if ((Diff_Mat[(*it1).first][(*it2).first] - LogE_unmerged) > Diff_LogE_best)
                {
                    Buffer_Community = MCM_unmerged;
                    Diff_LogE_best = Diff_Mat[(*it1).first][(*it2).first] - LogE_unmerged;
                    LogE_merged_best = Diff_Mat[(*it1).first][(*it2).first];
                }

                MCM_unmerged.clear();
            }
        }

        // ********* PERFORM MERGING *****:
        it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
        it++;  i_erase = (*it).first; //  new_part += (*it).second;

        for (it1 = Best_MCM_Partition.begin(); it1 != Best_MCM_Partition.end(); it1++)
        {
            Diff_Mat[(*it1).first][i_keep] = 0;
            Diff_Mat[i_keep][(*it1).first] = 0;
        }

        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; }
        iteration++;
    } // End While

    // ********* FINAL PARTITION **************:
    cout << endl << "******** FINAL PARTITION -- iteration " << iteration << ":" << endl << endl;

    LogE_tot = 0;
    for (it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
    {
        cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
        LogE_tot += Best_MCM_LogE_Part[(*it).first];
    }

    cout << endl << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
    cout << "\t LogE = " << LogE_tot << endl;   // LogE_MCM(Kset, Best_MCM_Partition, N) 

    return Best_MCM_Partition;
}

// Each iteration only considers 'k' communities to merge
map<unsigned int, uint64_t> Random(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    // Create initial independent model:
    map<unsigned int, uint64_t> Best_MCM_Partition, Best_MCM_copy, Merger_candidate;  uint64_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_Part;

    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = Op;
        Best_MCM_LogE_Part[i] = LogE_SubCM(Kset, Op, N);
        Op = Op << 1;
    }

    // Matrix to store change in LogE when communities are merged. (no double computations)
    double** Diff_Mat = (double**)malloc(r * sizeof(double*));

    for (unsigned int i = 0; i < r; i++)
    {
        Diff_Mat[i] = (double*)malloc(r * sizeof(double));
        for (int j = 0; j < r; j++)
        {
            Diff_Mat[i][j] = 0;
        }
    }

    // buffers:
    map<unsigned int, uint64_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, uint64_t>::iterator it1, it2, it;
    bool stop = false;

    int i_keep, i_erase; uint64_t new_part;   // define variables for merging

    int k = 1, iteration = 0;

    // ********* HIERARCHICAL MERGING ****************
    while (Best_MCM_Partition.size() > 1 && (!stop))
    {
        if (iteration > n / 2) { k = 2; }

        Diff_LogE_best = 0;
        LogE_merged = 0, LogE_unmerged = 0;

        // ***** Only consider two nodes to be merged with everything else ***** //
        Best_MCM_copy = Best_MCM_Partition;
        for (unsigned int i = 0; i < k; i++)
        {
            it = Best_MCM_copy.begin();
            advance(it, rand() % Best_MCM_copy.size());
            i_erase = it->first;
            Merger_candidate.insert(*it);
            Best_MCM_copy.erase(i_erase);
        }

        // ****** Find the best pair to merge: *********
        for (it1 = Merger_candidate.begin(); it1 != Merger_candidate.end(); it1++)
        {
            for (it2 = Best_MCM_Partition.begin(); it2 != Best_MCM_Partition.end(); it2++)
            {
                if ((*it1).first == (*it2).first) { continue; }
                //cout << counter << "\t" << bitset<n>((*it1).second) << "\t " << bitset<n>((*it2).second);

                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

                if (Diff_Mat[(*it1).first][(*it2).first] == 0)
                {
                    Diff_Mat[(*it1).first][(*it2).first] = LogE_SubCM(Kset, (*it1).second + (*it2).second, N) - LogE_unmerged;
                }
                //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged)

                if (Diff_Mat[(*it1).first][(*it2).first] > Diff_LogE_best)
                {
                    Buffer_Community = MCM_unmerged;
                    Diff_LogE_best = Diff_Mat[(*it1).first][(*it2).first];
                    LogE_merged_best = LogE_unmerged + Diff_Mat[(*it1).first][(*it2).first];
                }
                //cout << "\t " << Diff_LogE_best << endl;
                //counter++;
                MCM_unmerged.clear();
            }
        }

        for (it1 = Best_MCM_Partition.begin(); it1 != Best_MCM_Partition.end(); it1++)
        {
            for (it2 = Best_MCM_Partition.begin(); it2 != Best_MCM_Partition.end(); it2++)
            {
                MCM_unmerged.insert(*it1);
                MCM_unmerged.insert(*it2);
                if (Diff_Mat[(*it1).first][(*it2).first] > Diff_LogE_best)
                {
                    Diff_LogE_best = Diff_Mat[(*it1).first][(*it2).first];
                    Buffer_Community = MCM_unmerged;
                }
                MCM_unmerged.clear();
            }
        }

        // ********* PERFORM MERGING *****:
        it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
        it++;  i_erase = (*it).first; //  new_part += (*it).second;

        for (it1 = Best_MCM_Partition.begin(); it1 != Best_MCM_Partition.end(); it1++)
        {
            Diff_Mat[(*it1).first][i_keep] = 0;
            Diff_Mat[i_keep][(*it1).first] = 0;
        }

        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; }

        Merger_candidate.clear();
        iteration++;
    } // End While

    // ********* FINAL PARTITION **************:
    cout << endl << "******** FINAL PARTITION -- iteration " << iteration << ":" << endl << endl;

    LogE_tot = 0;
    for (it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
    {
        cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
        LogE_tot += Best_MCM_LogE_Part[(*it).first];
    }

    cout << endl << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
    cout << "\t LogE = " << LogE_tot << endl;   // LogE_MCM(Kset, Best_MCM_Partition, N) 

    return Best_MCM_Partition;
}

// Use maximum from matrix and exclude random communities that have already been computed.
map<unsigned int, uint64_t> Random2(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    map<unsigned int, uint64_t> Best_MCM_Partition, MCM_part;  uint64_t Op = 1;
    map<unsigned int, double> Best_MCM_LogE_Part;

    map<unsigned int, map<unsigned int, double>> Diff_mat;
    vector<double> vec;
    vector<unsigned int> candidate;

    // Initial independent model
    for (unsigned int i = 0; i < r; i++)
    {
        Best_MCM_Partition[i] = Op;
        Best_MCM_LogE_Part[i] = LogE_SubCM(Kset, Op, N);
        Op = Op << 1;

        // Matrix of mergers
        for (unsigned int j = i+1; j < r; j++)
        {
            Diff_mat[j][i] = 0.0;
        }

        // k communities to choose from
        candidate.push_back(i);
    }

    // iterators
    map<unsigned int, uint64_t>::iterator it, it_part;
    vector<unsigned int>::iterator it_vec;
    // Log evidence values
    double Diff_LogE_best, LogE_unmerged, LogE_merged, LogE_merged_best;
    // indices, number of random communities chosen, iteration count, stopping variable
    int end, i_keep, i_erase, n1, n2, k = 1, iteration = 1;
    bool stop = false, erased = false, kept = false;

    int control;

    while (Best_MCM_Partition.size() > 1 && (!stop))
    {

        // tuning of k
        if (iteration >= r / 2) { k = 2; }

        Diff_LogE_best = 0;

        // ***** Shuffle vector and use first k values to get random communities ***** //
        random_shuffle(candidate.begin(), candidate.end());
        if ( ((signed)candidate.size() - k) >= 0) 
        { 
            end = ((signed)candidate.size() - k);
        } 
        else
        { 
            end = 0; 
        }
        
        for (unsigned int i = candidate.size(); i > end; i--)
        {
            MCM_part[candidate[i-1]] = Best_MCM_Partition[candidate[i-1]];
            it_part = MCM_part.begin();
            for (it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
            {
                if ((*it_part).first == (*it).first) { continue; }

                if ((*it_part).first > (*it).first) 
                { 
                    n1 = (*it_part).first; 
                    n2 = (*it).first; 
                }
                else 
                { 
                    n1 = (*it).first; 
                    n2 = (*it_part).first; 
                }

                LogE_unmerged = Best_MCM_LogE_Part[n1] + Best_MCM_LogE_Part[n2];

                if (Diff_mat[n1][n2] == 0)
                {
                    LogE_merged = LogE_SubCM(Kset, (*it_part).second + (*it).second, N);
                    Diff_mat[n1][n2] = LogE_merged - LogE_unmerged;
                }

                if (Diff_mat[n1][n2] > Diff_LogE_best)
                {
                    Diff_LogE_best = Diff_mat[n1][n2];
                }
            }
            MCM_part.clear();
        }
        //cout << "Iteration: " << iteration << endl;
        //cout << "Candidate size: " << candidate.size() << endl;
        //print_values(Diff_mat);
        //print_keys(Diff_mat);

        // ********* GET MAX AND INDEX **************:
        Diff_LogE_best = return_max(Diff_mat, i_keep, i_erase);
        LogE_merged_best = Best_MCM_LogE_Part[i_keep] + Best_MCM_LogE_Part[i_erase] + Diff_LogE_best;

        //cout << "keep: " << i_keep << ". erase: " << i_erase << "." << endl << endl;

        for (auto& pair : Diff_mat[i_keep]) 
        { 
            pair.second = 0; // Set the values in the 'row' of i_keep to 0
        }

        for (auto& pair : Diff_mat) 
        { 
            Diff_mat[pair.first].erase(i_erase); // Erase the 'column' of i_erase
            if (Diff_mat[pair.first].count(i_keep)) 
            { 
                Diff_mat[pair.first][i_keep] = 0; // Set the values in the 'column' of i_erase to 0
            } 
        }
        
        Diff_mat.erase(i_erase); // Erase the 'row' of i_erase

        // ********* PERFORM MERGING ****************:
        Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
        Best_MCM_Partition.erase(i_erase);

        Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
        Best_MCM_LogE_Part.erase(i_erase);

        // ********* EXCLUDE CANDIDATES *************:
        if (k < (signed)candidate.size()) { end = k; }
        else { end = (signed)candidate.size(); }

        for (unsigned int i = 0; i < end; i++)
        {
            if (candidate.back() == i_keep) { kept = true; }
            if (candidate.back() == i_erase) { erased = true; }
            candidate.pop_back();
        }

        if (kept || candidate.empty()) { candidate.push_back(i_keep); kept = false; }

        if (!erased && candidate.size()!= 1)
        {
            for (it_vec = candidate.begin(); it_vec != candidate.end(); it_vec++)
            {
                if (*it_vec == i_erase) { candidate.erase(it_vec); break; }
            }
        }
        erased = false;

        // ********* STOPPING CRITERIA **************:
        if (Diff_LogE_best == 0) { stop = true; }

        /*
        cout << endl << "Enter '1' if you want to exit: ";
        cin >> control;
        if (control == 1) { break; }
        */
        iteration++;
    }
    return Best_MCM_Partition;
}

// Variation of information
double Var_of_Inf( map<unsigned int, uint64_t> Partition1, map<unsigned int, uint64_t> Partition2 )
{
    // Variation of information calculates the distance between two partitions. The regular variation of information
    // is equal to the joint entropy minus the mutual information. However, the normalized version (divide by the joint entrop)
    // is preferred over the regular as this is a true metric, i.e., it satisfies the triangle inequality.
    double I, H, p1, p2, p12;
    I = 0;
    H = 0;
    map<unsigned int, uint64_t>::iterator com1, com2;
    for (com1 = Partition1.begin(); com1 != Partition1.end(); com1++)
    {
        p1 = (double)(bitset<n>((*com1).second).count()) / (double)(n);
        for (com2 = Partition2.begin(); com2 != Partition2.end(); com2++)
        {
            p2 = (double)(bitset<n>((*com2).second).count()) / (double)(n);
            p12 = (double)(bitset<n>((*com1).second & (*com2).second).count()) / (double)(n);
            if (p12 != 0)
            {
                I += p12 * log(p1 * p2);
                H += p12 * log(p12);
            }
        }
    }
    if(H == 0) { return 0; }
    else { return 2 - I / H; }
}

/******************************************************************************/
/*******************************   main function   ****************************/
/******************************************************************************/
int main()
{
    cout << "--->> Create OUTPUT Folder: (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "******************************  Choice of the basis:  *************************************";
    cout << endl << "*******************************************************************************************" << endl;
    // *** Choice of the basis for building the Minimally Complex Model (MCM):

    // *** Basis elements are written using the integer representation of the operator
    // *** For instance, a basis element on the last two spin variable would be written: 
    // ***      -->  Op = s1 s2           Spin operator
    // ***      -->  Op = 000000011       Binary representation
    // ***      -->  Op = 3               Integer representation   ( 000000011 = 3 )

      // *** The basis can be specified by hand here:
      /*uint64_t Basis_Choice[] = {36, 10, 3, 272, 260, 320, 130, 65, 4};    // Ex. This is the best basis for the USSC dataset

      unsigned int m = sizeof(Basis_Choice) / sizeof(uint64_t);
      list<uint64_t> Basis_li;  Basis_li.assign (Basis_Choice, Basis_Choice + m);
      */

      // *** Basis can also be read from a file:
      // list<uint64_t> Basis_li = Read_BasisOp_IntegerRepresentation();

      // *** Or one can use the original basis:
    list<uint64_t> Basis_li = Original_Basis();

    // *** Print info about the Basis:
    PrintTerm_Basis(Basis_li);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  Read the data:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;
    unsigned int N = 0; // will contain the number of datapoints in the dataset
    map<uint64_t, unsigned int> Nset = read_datafile(&N);

    cout << endl << "*******************************************************************************************";
    cout << endl << "*********************************  Change the data basis   ********************************";
    cout << endl << "**************************************  Build Kset:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;
    // *** Transform the data in the specified in Basis_SCModel[];
    map<uint64_t, unsigned int> Kset = build_Kset(Nset, Basis_li, false);

    cout << endl << "*******************************************************************************************";
    cout << endl << "************************************  All Indep Models:  **********************************";
    cout << endl << "*******************************************************************************************" << endl;

    cout << "Independent models in the new basis:" << endl;

    unsigned int r = 20;

    map<unsigned int, uint64_t> Partition_Indep;  uint64_t Op = 1;
    for (unsigned int i = 0; i < r; i++)
    {
        Partition_Indep[i] = Op;
        //cout << "Added Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;    
        Op = Op << 1;
    }
    cout << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;

    /*
    cout << endl << "*******************************************************************************************";
    cout << endl << "****************************************  Original:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;
    // ************ Original:
    auto start = chrono::system_clock::now();
    Original(Kset, N);

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed1 = end - start;
    cout << endl << "Elapsed time: " << elapsed1.count() << "s";
    */
    
    

    cout << endl << "*******************************************************************************************";
    cout << endl << "*****************************************  Matrix:  ***************************************";
    cout << endl << "*******************************************************************************************" << endl;
    // ************ Matrix:
    map<unsigned int, uint64_t> MCM_Partition;

    auto start = chrono::system_clock::now();
    MCM_Partition = Matrix(Kset, N);

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed2 = end - start;
    cout << endl << "Elapsed time: " << elapsed2.count() << "s" << endl << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "*****************************************  Random:  ***************************************";
    cout << endl << "*******************************************************************************************" << endl;
    // ************ Random:
    srand(time(NULL));
    map<unsigned int, uint64_t> Random_Partition;
    map<unsigned int, uint64_t>::iterator iter;
    vector<double> VI, TIME;
    vector<double>::iterator it;

    double voi, z = 1.96;
    int nSamples = 50;
    
    /*for (iter = Random_Partition.begin(); iter != Random_Partition.end(); iter++)
    {
        cout << (*iter).first << "\t " << bitset<n>((*iter).second) << endl;
    }*/

    
    for (unsigned int i = 0; i < nSamples; i++)
    {
        start = chrono::system_clock::now();
        Random_Partition = Random2(Kset, N);
        end = chrono::system_clock::now();

        elapsed2 = end - start;
        cout << elapsed2.count() << ", ";

        voi = Var_of_Inf(MCM_Partition, Random_Partition);
        VI.push_back(voi);
    }


    cout << endl << endl;

    double sum = accumulate(VI.begin(), VI.end(), 0.0);
    double mean = sum / VI.size();
    sum = inner_product(VI.begin(), VI.end(), VI.begin(), 0.0);
    double std = sqrt(sum / VI.size() - mean * mean);

    cout << endl << "Mean: " << mean << endl << "Std: " << std;
    cout << endl << "Conf: [" << mean - z * std / (double)nSamples << ", " << mean + z * std / (double)nSamples << ']' << endl;

    return 0;
}