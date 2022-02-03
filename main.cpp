//g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Basis_Choice.cpp RandomPart.cpp 
// Best_MCM.cpp 
//time ./a.out
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

//map<unsigned int, uint64_t> Random_Partition(unsigned int A, mt19937 gen);

unsigned int A = 5;   // A = number of parts
uniform_int_distribution<int> uni(0., A - 1);

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

map<unsigned int, uint64_t> Matrix(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
    // Create initial independent model:
    map<unsigned int, uint64_t> Best_MCM_Partition;  uint64_t Op = 1;
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
    map<unsigned int, uint64_t> MCM_unmerged, Buffer_Community;
    double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0, LogE_merged_best = 0, LogE_tot = 0;

    // Best:
    map<unsigned int, uint64_t>::iterator it1, it2, it1_end, it2_start, it;
    int counter = 0;
    bool stop = false;

    int i_keep, i_erase, i_mat; uint64_t new_part;   // define variables for merging

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


double Var_of_Inf(map<unsigned int, uint64_t> Partition1, map<unsigned int, uint64_t> Partition2)
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
    if (H == 0) { return 0; }
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

    // *** Specify path to data in 'data.h'.

    unsigned int N = 0; // will contain the number of datapoints in the dataset
    map<uint64_t, unsigned int> Nset = read_datafile(&N);

    if (N == 0) { return 0; } // Terminate program if the file can't be found

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

    unsigned int r = n;

    map<unsigned int, uint64_t> Partition_Indep;  uint64_t Op = 1;
    for (unsigned int i = 0; i < r; i++)
    {
        Partition_Indep[i] = Op;
        cout << "Added Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;    
        Op = Op << 1;
    }
    cout << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "******************************  Hierachical merging result:  ******************************";
    cout << endl << "*******************************************************************************************" << endl;

    map<unsigned int, uint64_t> Final_Partition;
    // *** Calculate the optimal partition
    auto start = chrono::system_clock::now();
    Final_Partition = Matrix(Kset, N);
    auto end = chrono::system_clock::now();

    // *** Check if the partition is valid
    if (!check_partition(Final_Partition)) { cout << "Result is invalid! A node appears in more than one community!" << endl << endl; }

    // *** Time it takes to find partition
    chrono::duration<double> elapsed = end - start;
    cout << endl << "Elapsed time: " << elapsed.count() << "s" << endl << endl;

    return 0;
}
