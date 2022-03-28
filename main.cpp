// To compile: g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Basis_Choice.cpp RandomPart.cpp 
// To run: time ./a.out
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <bitset>
#include <map>
#include <cmath>       /* tgamma */
#include <random>

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
std::mt19937 gen;       //only for mt19937    //std::default_random_engine gen;
void initialise_generator()
{
    int seed = (unsigned)time(NULL);
    //srand48(seed);      //for drand48
    gen.seed(seed);     //for mt19937
}

unsigned int A = 5;   // A = number of parts
uniform_int_distribution<int> uni(0., A - 1);

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

map<unsigned int, __int128_t> Matrix(map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
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

    int i_keep, i_erase, i_mat; __int128_t new_part;   // define variables for merging

    int k = 0, iteration = 0;

    // ********* HIERARCHICAL MERGING ****************
    while (Best_MCM_Partition.size() > 1 && (!stop))
    {
        //prints: 
        /*= Best_MCM_Partition.size();
        cout << endl << "******** Current partition -- iteration " << iteration << ":" << endl;
        LogE_tot = 0;
        for(it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
        {
          cout << "\t" << (*it).first << "\t " << bitset<n>((*it).second) << endl;
          LogE_tot += Best_MCM_LogE_Part[(*it).first];
        }
        cout << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
        cout << "\t LogE = " << LogE_tot << endl;     // LogE_MCM(Kset, Best_MCM_Partition, N)
        */
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

        // ********* PRINTS **************:
        //cout << endl << "Counter = " << counter << "\t VS \t nb of pairs = " << k*(k-1)/2 << endl << endl;
        //cout << "Best Merging: LogE = " << Diff_LogE_best << endl;
        //cout << "Nb of pairs = " << k * (k - 1) / 2 << endl;
        /*
        for(it = Buffer_Community.begin(); it != Buffer_Community.end(); it++)
        {
          cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
        } */

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

/******************************************************************************/
/*******************************   main function   ****************************/
/******************************************************************************/
int main()
{
    cout << "--->> Create OUTPUT Folder: (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************************  Sample data:  *****************************************";
    cout << endl << "*******************************************************************************************" << endl;

    //system(("C:/Users/sam/source/repos/MinCompSpin_Greedy_Private-main/INPUT/binary_networks/benchmark -f C:/Users/sam/source/repos/MinCompSpin_Greedy_Private-main/INPUT/binary_networks/flags.dat -mu " + to_string(0.2) + " > nul:").c_str());
    list<Interaction> list_I = write_interactions(0.2, "datasets/network113.dat");
    Sample_dataset(list_I, datafilename, 10000);

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
      /*__int128_t Basis_Choice[] = {36, 10, 3, 272, 260, 320, 130, 65, 4};    // Ex. This is the best basis for the USSC dataset

      unsigned int m = sizeof(Basis_Choice) / sizeof(__int128_t);
      list<__int128_t> Basis_li;  Basis_li.assign (Basis_Choice, Basis_Choice + m);
      */

      // *** Basis can also be read from a file:
      // list<__int128_t> Basis_li = Read_BasisOp_IntegerRepresentation();

      // *** Or one can use the original basis:
    list<__int128_t> Basis_li = Original_Basis();

    // *** Print info about the Basis:
    //PrintTerm_Basis(Basis_li);

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  Read the data:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;

    unsigned int N = 0; // will contain the number of datapoints in the dataset
    map<__int128_t, unsigned int> Kset = read_datafile(&N);

    if (N == 0) { return 0; } // Terminate program if the file can't be found

    //cout << endl << "*******************************************************************************************";
    //cout << endl << "*********************************  Change the data basis   ********************************";
    //cout << endl << "**************************************  Build Kset:  **************************************";
    //cout << endl << "*******************************************************************************************" << endl;
    //// *** Transform the data in the specified in Basis_SCModel[];
    ////map<__int128_t, unsigned int> Kset = build_Kset(Nset, Basis_li, false);

    cout << lgamma((double)(N + (1UL << (10 - 1))));

    cout << endl << "*******************************************************************************************";
    cout << endl << "******************************  Hierachical merging result:  ******************************";
    cout << endl << "*******************************************************************************************" << endl;

    map<unsigned int, __int128_t> fp1, fp2;
    // *** Calculate the optimal partition
    auto start = chrono::system_clock::now();
    fp1 = Matrix(Kset, N);
    auto end = chrono::system_clock::now();

    // *** Check if the partition is valid
    if (!check_partition(fp1)) { cout << "Result is invalid! A node appears in more than one community!" << endl << endl; }

    // *** Time it takes to find partition
    chrono::duration<double> elapsed = end - start;
    
    cout << "######### EMPERICAL #########" << endl;
    double H = Entropy(Kset, N);
    cout << "H : " << H << ". Range: [0, " << N << "]" << endl << endl;

    cout << "#########  GREEDY   #########" << endl;
    double LE = LogE_MCM(Kset, fp1, N);
    Print_Partition(fp1);
    
    map<__int128_t, double> greedy_prob = MCM_distr(Kset, fp1, N);
    map<__int128_t, double> emp_prob = emp_dist(Kset, N);
    double JSD = JS_divergence(greedy_prob, emp_prob, N);

    cout << "Elapsed time    : " << elapsed.count() << "s" << endl;
    cout << "Log-evidence    : " << LE << endl;
    cout << "JSD(emp || mcm) : " << JSD << endl << endl;

    //**********************************************
    //******** MATRIX 2 ****************************
    //**********************************************

    //start = chrono::system_clock::now();
    //fp2 = Matrix2(Kset, N);
    //end = chrono::system_clock::now();

    //// *** Check if the partition is valid
    //if (!check_partition(fp2)) { cout << "Result is invalid! A node appears in more than one community!" << endl << endl; }

    //// *** Time it takes to find partition
    //elapsed = end - start;
    //cout << endl << "Elapsed time: " << elapsed.count() << "s" << endl << endl;

    fp2 = read_communities("datasets/community113.dat");

    cout << "#########   EXACT   #########" << endl;
    Print_Partition(fp2);

    map<__int128_t, double> exact_prob = MCM_distr(Kset, fp2, N);
    LE = LogE_MCM(Kset, fp2, N);
    cout << "Log-evidence    : " << LE << endl;
    
    JSD = JS_divergence(exact_prob, emp_prob, N);
    cout << "JSD(emp || mcm) : " << JSD << endl << endl;

    JSD = JS_divergence(exact_prob, greedy_prob, N);

    cout << "####### COMPARE GREEDY OUTCOUME TO EXACT RESULT #######" << endl;
    cout << "Jensen-Shannon divergence:     " << JSD << endl;

    cout << "Variation of information:      " << Var_of_Inf(fp1, fp2) << endl;
    cout << "Normalized mutual information: " << Norm_Mut_info(fp1, fp2) << endl << endl;

    return 0;
}