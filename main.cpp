// To compile: g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp info_quant.cpp best_basis.cpp Basis_Choice.cpp P_s.cpp metropolis.cpp Best_MCM_Greedy.cpp
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
using namespace std::chrono;

/******************************************************************************/
/**********************    CONSTANTS AND FUNCTIONS    *************************/
/******************************************************************************/
#include "data.h"
#include "library.h"
//#include "library_Metropolis.h"

/******************************************************************************/
/*******************************   main function   ****************************/
/******************************************************************************/
int main()
{
    cout << "--->> Create OUTPUT Folder: (if needed) ";
    system(("mkdir -p " + OUTPUT_directory).c_str());
    cout << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "***********************************  READ THE DATA:  **************************************";
    cout << endl << "*******************************************************************************************" << endl;

    unsigned int N = 0; // will contain the number of datapoints in the dataset
    map<__int128_t, unsigned int> Nset = read_datafile(&N);

    if (N == 0) { return 0; } // Terminate program if the file can't be found

    cout << endl << "                        ###### File has been read succesfully ######" << endl;
    cout << "Number of datapoints = " << N << endl;
    cout << "Number of different observed states = " << Nset.size() << endl;


    cout << endl << "*******************************************************************************************";  
    cout << endl << "******************************  CHOICE OF THE BASIS:  *************************************";
    cout << endl << "*******************************************************************************************" << endl;

//    list<__int128_t> Basis_li = Original_Basis();  // original basis of the data: this is the most natural choice a priori

  // *** The basis can also be read from a file:
    list<__int128_t> Basis_li = Read_BasisOp_IntegerRepresentation();
//   list<__int128_t> Basis_li = Read_BasisOp_BinaryRepresentation();

    PrintTerm_Basis(Basis_li);


    cout << endl << "*******************************************************************************************";
    cout << endl << "************************  Transform the data in the new basis   ***************************";
    cout << endl << "**********************************   Build Kset:   ****************************************";
    cout << endl << "*******************************************************************************************" << endl;
    //// *** Transform the data in the specified in Basis_SCModel[];
    map<__int128_t, unsigned int> Kset = build_Kset(Nset, Basis_li, false);


    cout << endl << "*******************************************************************************************";
    cout << endl << "******************************  Hierachical merging result:  ******************************";
    cout << endl << "*******************************************************************************************" << endl;

    // *** Calculate the optimal partition
    auto start = chrono::system_clock::now();
    map<unsigned int, __int128_t> fp1 = Matrix(Kset, N);
    auto end = chrono::system_clock::now();

    // *** Time it takes to find partition
    chrono::duration<double> elapsed = end - start;

    cout << "######### EMPERICAL #########" << endl;
    // Entropy of dataset
    double H = Entropy(Kset, N);
    cout << "H : " << H << ". Range: [0, " << n << "]" << endl << endl;

    cout << "#########  GREEDY   #########" << endl;
    // Log evidence of MCM
    double LE_g = LogE_MCM(Kset, fp1, N);
    Print_Partition(fp1);

    cout << "Elapsed time      : " << elapsed.count() << "s" << endl;
    cout << "Log-evidence      : " << LE_g << endl;
    cout << "Average comm size : " << (double)n / (double)fp1.size() << endl << endl;

    cout << "#########  THEORETICAL   #########" << endl;
    map<unsigned int, __int128_t> fp2 = read_communities(communityfile);

    double LE_t = LogE_MCM(Kset, fp2, N);
    Print_Partition(fp2);

    cout << "Log-evidence      : " << LE_t << endl;
    cout << "Average comm size : " << (double)n / (double)fp2.size() << endl << endl;

    cout << "#########  COMPARATIVE MEASURES   #########" << endl;
    double VOI = Var_of_Inf(fp1, fp2);
    double NMI = Norm_Mut_info(fp1, fp2);
    string istrue = is_subset(fp1, fp2) ? "Yes" : "No";

    cout << "Is MCM_g \'subset\' of MCM_t    : " << istrue << endl;
    cout << "Variation of Information      : " << VOI << endl;
    cout << "Normalized Mutual Information : " << NMI << endl;
    cout << "Difference in Log-Evidence    : " << LE_g - LE_t << endl << endl;

    cout << endl << "*******************************************************************************************";
    cout << endl << "**********************  Print information about the found MCM:  ***************************";
    cout << endl << "*******************************************************************************************" << endl;

    // Prints 1) information about the MCM; 2) the state probabilities P(s) (in the Data VS MCM); 3) the probability P(k) of observing a state with k values "+1" (in the Data VS MCM) 
    PrintFile_StateProbabilites_OriginalBasis(Nset, Basis_li, fp1, N, "Result");

    // Print the state probabilities P(s) (in the Data VS MCM) using the data transformed in the bew basis:
    PrintFile_StateProbabilites_NewBasis(Kset, fp1, N, "Result");

    return 0;
}