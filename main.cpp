//g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Basis_Choice.cpp RandomPart.cpp //Best_MCM.cpp 
//time ./a.out
//
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
          srand48(seed);      //for drand48
          gen.seed(seed);     //for mt19937
}

//map<unsigned int, uint64_t> Random_Partition(unsigned int A, mt19937 gen);

unsigned int A=5;   // A = number of parts
uniform_int_distribution<int> uni(0.,A-1);

map<unsigned int, uint64_t> Random_Partition()
{
//  uniform_int_distribution<int> uni(0.,A-1);
  map<unsigned int, uint64_t> Partition_rand;
  uint64_t si = 0;
  unsigned int part = 0;

  for(unsigned int a=0; a<A; a++)
  {
    Partition_rand[a]=si;
  }

  si = 1;

  for(int i=0; i<n; i++)
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
map<uint32_t, uint32_t> Convert_Partition_forMCM(uint32_t *a, unsigned int r = n)
{
  map<uint32_t, uint32_t> Partition;
  uint32_t element = 1;

  for (int i=r-1; i>=0; i--)  // read element from last to first
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

double LogE_MCM_converted(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int *a, unsigned int r = n)
{
  map<unsigned int, uint64_t> Partition;
  uint64_t element = 1;

  for (int i=r-1; i>=0; i--)  // read element from last to first
    {    
      Partition[(a[i])] += element;
      element = element << 1;      //cout << a[i] << "\t :";   Print_Partition_Converted(Partition);
    }

  double LogE = LogE_MCM(Kset, Partition, N);
  Partition.clear();

  return LogE;
}

void Print_Partition(uint32_t *a, unsigned int r = n)
{
  //for (int i=0; i<n-r; i++)
  //{   cout << "x" << " ";   }
  for (int i=0; i<r; i++)
  {    cout << a[i] << " ";  }
  cout << endl;
}

// map<unsigned int, uint64_t>
void Hierarchical_merging(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
  // Create initial independent model:
  unsigned int *a_best = (unsigned int *)malloc(r*sizeof(unsigned int));

  for (unsigned int i=0; i<r; i++)
  {    a_best[i]=i;  }

  // LogE of Best_MCM:
  Print_Partition(a_best, r);
  double LogE_best = LogE_MCM_converted(Kset, N, a_best, r);
  cout << endl;


  cout << LogE_best << endl;
}

// void PrintTerminal_MCM_Info(map<uint64_t, unsigned int> Kset, unsigned int N, map<unsigned int, uint64_t> MCM_Partition);

//map<unsigned int, uint64_t> 

void Hierarchical_merging2(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
  // Create initial independent model:
  map<unsigned int, uint64_t> Best_MCM_Partition;  uint64_t Op = 1;
  for (unsigned int i = 0 ; i<r; i++)
  {  
    Best_MCM_Partition[i] = Op;
    Op = Op << 1;
  }

  // LogE of Best_MCM_Partition:
  double Best_LogE = LogE_MCM(Kset, Best_MCM_Partition, N);
  //cout << Best_LogE << endl;

  // buffers:
  map<unsigned int, uint64_t> MCM_merged, MCM_unmerged, Buffer_Community;
  double Diff_LogE_best = 0, LogE_merged = 0, LogE_unmerged = 0;

  // Best:
  map<unsigned int, uint64_t>::iterator it1, it2, it1_end, it2_start, it;
  int counter = 0;
  bool stop = false;

  int i_keep, i_erase; uint64_t new_part;   // define variables for merging

  int k=0, iteration = 0;

  // ********* HIERARCHICAL MERGING ****************
  while(Best_MCM_Partition.size() > 1 && (!stop) )
  { 
    //prints: 
    k= Best_MCM_Partition.size();
    cout << endl << "******** Current partition -- iteration " << iteration << ":" << endl;
    for(it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
    {
      cout << "\t" << (*it).first << "\t " << bitset<n>((*it).second) << endl;
    } 
    cout << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
    cout << "\t LogE = " << LogE_MCM(Kset, Best_MCM_Partition, N) << endl;

    counter = 0;
    Diff_LogE_best = 0;
    LogE_merged = 0, LogE_unmerged = 0;

    // ****** Find the best pair to merge: *********
    //cout << endl << "******** Test combinations: " << endl;
    it1_end=Best_MCM_Partition.end();   it1_end--;
    for(it1 = Best_MCM_Partition.begin(); it1 != it1_end; it1++)
    {    
      it2_start = it1;    it2_start++;
      for(it2 = it2_start; it2!= Best_MCM_Partition.end(); it2++)
      {
        //cout << counter << "\t" << bitset<n>((*it1).second) << "\t " << bitset<n>((*it2).second);
      
        MCM_unmerged.insert(*it1);
        MCM_unmerged.insert(*it2);
        LogE_unmerged = LogE_MCM(Kset, MCM_unmerged, N);

        MCM_merged[0]=(*it1).second + (*it2).second;
        LogE_merged = LogE_MCM(Kset, MCM_merged, N);

        //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged);

        if( (LogE_merged-LogE_unmerged) > Diff_LogE_best)
        {
          Buffer_Community = MCM_unmerged;
          Diff_LogE_best = LogE_merged-LogE_unmerged;
        }
        //cout << "\t " << Diff_LogE_best << endl;
        //counter++;
        MCM_unmerged.clear();
      }
    }

    // ********* PRINTS **************:
    //cout << endl << "Counter = " << counter << "\t VS \t nb of pairs = " << k*(k-1)/2 << endl << endl;
    cout << "Best Merging: LogE = " << Diff_LogE_best << endl;
    cout << "Nb of pairs = " << k*(k-1)/2 << endl;
    /*
    for(it = Buffer_Community.begin(); it != Buffer_Community.end(); it++)
    {
      cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
    } */

    // ********* PERFORM MERGING *****:
    it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
    it++;  i_erase = (*it).first; //  new_part += (*it).second;

    Best_MCM_Partition[i_keep]=Best_MCM_Partition[i_keep]+Best_MCM_Partition[i_erase];
    Best_MCM_Partition.erase(i_erase);

    // ********* STOPPING CRITERIA **************:
    if (Diff_LogE_best == 0) {  stop = true; }
    iteration++;
  } // End While

  // ********* FINAL PARTITION **************:
  cout << endl << "******** FINAL PARTITION -- iteration " << iteration << ":" << endl << endl;

  for(it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
  {
    cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
  } 

  cout << endl << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
  cout << "\t LogE = " << LogE_MCM(Kset, Best_MCM_Partition, N) << endl;

}


void Hierarchical_merging3(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
  // Create initial independent model:
  map<unsigned int, uint64_t> Best_MCM_Partition;  uint64_t Op = 1;
  map<unsigned int, double> Best_MCM_LogE_Part;

  for (unsigned int i = 0 ; i<r; i++)
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

  int k=0, iteration = 0;

  // ********* HIERARCHICAL MERGING ****************
  while(Best_MCM_Partition.size() > 1 && (!stop) )
  { 
    //prints: 
    k= Best_MCM_Partition.size();
    cout << endl << "******** Current partition -- iteration " << iteration << ":" << endl;
    LogE_tot = 0;
    for(it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
    {
      cout << "\t" << (*it).first << "\t " << bitset<n>((*it).second) << endl;
      LogE_tot += Best_MCM_LogE_Part[(*it).first];
    } 
    cout << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
    cout << "\t LogE = " << LogE_tot << endl;     // LogE_MCM(Kset, Best_MCM_Partition, N) 

    counter = 0;
    Diff_LogE_best = 0;
    LogE_merged = 0, LogE_unmerged = 0;

    // ****** Find the best pair to merge: *********
    //cout << endl << "******** Test combinations: " << endl;
    it1_end=Best_MCM_Partition.end();   it1_end--;
    for(it1 = Best_MCM_Partition.begin(); it1 != it1_end; it1++)
    {    
      it2_start = it1;    it2_start++;
      for(it2 = it2_start; it2!= Best_MCM_Partition.end(); it2++)
      {
        //cout << counter << "\t" << bitset<n>((*it1).second) << "\t " << bitset<n>((*it2).second);
      
        MCM_unmerged.insert(*it1);
        MCM_unmerged.insert(*it2);
        LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

        LogE_merged = LogE_SubCM(Kset, (*it1).second + (*it2).second, N);

        //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged);

        if( (LogE_merged-LogE_unmerged) > Diff_LogE_best)
        {
          Buffer_Community = MCM_unmerged;
          Diff_LogE_best = LogE_merged-LogE_unmerged;
          LogE_merged_best = LogE_merged;
        }
        //cout << "\t " << Diff_LogE_best << endl;
        //counter++;
        MCM_unmerged.clear();
      }
    }

    // ********* PRINTS **************:
    //cout << endl << "Counter = " << counter << "\t VS \t nb of pairs = " << k*(k-1)/2 << endl << endl;
    cout << "Best Merging: LogE = " << Diff_LogE_best << endl;
    cout << "Nb of pairs = " << k*(k-1)/2 << endl;
    /*
    for(it = Buffer_Community.begin(); it != Buffer_Community.end(); it++)
    {
      cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
    } */

    // ********* PERFORM MERGING *****:
    it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
    it++;  i_erase = (*it).first; //  new_part += (*it).second;

    Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
    Best_MCM_Partition.erase(i_erase);

    Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
    Best_MCM_LogE_Part.erase(i_erase);

    // ********* STOPPING CRITERIA **************:
    if (Diff_LogE_best == 0) {  stop = true; }
    iteration++;
  } // End While

  // ********* FINAL PARTITION **************:
  cout << endl << "******** FINAL PARTITION -- iteration " << iteration << ":" << endl << endl;

  LogE_tot = 0;
  for(it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
  {
    cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
    LogE_tot += Best_MCM_LogE_Part[(*it).first];
  } 

  cout << endl << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
  cout << "\t LogE = " << LogE_tot << endl;   // LogE_MCM(Kset, Best_MCM_Partition, N) 
}


void Hierarchical_merging4(map<uint64_t, unsigned int> Kset, unsigned int N, unsigned int r = n)
{
  // Create initial independent model:
  map<unsigned int, uint64_t> Best_MCM_Partition;  uint64_t Op = 1;
  map<unsigned int, double> Best_MCM_LogE_Part;

  double** Diff_logE_Mat = (double**)malloc(r*sizeof(double *));

  for(int i=0; i<r; i++)
  {   
    Diff_logE_Mat[i] = (double*)malloc(r*sizeof(double));
    for(int j=0; j<r; j++)
      {   Diff_logE_Mat[i][j] = 0;  }
  }

  for (unsigned int i = 0 ; i<r; i++)
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

  int k=0, iteration = 0;

  // ********* HIERARCHICAL MERGING ****************
  while(Best_MCM_Partition.size() > 1 && (!stop) )
  { 
    //prints: 
    k= Best_MCM_Partition.size();
    cout << endl << "******** Current partition -- iteration " << iteration << ":" << endl;
    LogE_tot = 0;
    for(it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
    {
      cout << "\t" << (*it).first << "\t " << bitset<n>((*it).second) << endl;
      LogE_tot += Best_MCM_LogE_Part[(*it).first];
    } 
    cout << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
    cout << "\t LogE = " << LogE_tot << endl;     // LogE_MCM(Kset, Best_MCM_Partition, N) 

    counter = 0;
    Diff_LogE_best = 0;
    LogE_merged = 0, LogE_unmerged = 0;

    // ****** Find the best pair to merge: *********
    //cout << endl << "******** Test combinations: " << endl;
    it1_end=Best_MCM_Partition.end();   it1_end--;
    for(it1 = Best_MCM_Partition.begin(); it1 != it1_end; it1++)
    {    
      it2_start = it1;    it2_start++;
      for(it2 = it2_start; it2!= Best_MCM_Partition.end(); it2++)
      {
        //cout << counter << "\t" << bitset<n>((*it1).second) << "\t " << bitset<n>((*it2).second);
      
        MCM_unmerged.insert(*it1);
        MCM_unmerged.insert(*it2);
        LogE_unmerged = Best_MCM_LogE_Part[(*it1).first] + Best_MCM_LogE_Part[(*it2).first]; //LogE_MCM(Kset, MCM_unmerged, N);

        LogE_merged = LogE_SubCM(Kset, (*it1).second + (*it2).second, N);

        //cout << "\t DlogE = " << (LogE_merged-LogE_unmerged);
        Diff_logE_Mat[(*it1).first][(*it2).first] = ((LogE_merged-LogE_unmerged)>0)?(LogE_merged-LogE_unmerged):0;

        if( (LogE_merged-LogE_unmerged) > Diff_LogE_best)
        {
          Buffer_Community = MCM_unmerged;
          Diff_LogE_best = LogE_merged-LogE_unmerged;
          LogE_merged_best = LogE_merged;
        }
        //cout << "\t " << Diff_LogE_best << endl;
        //counter++;
        MCM_unmerged.clear();
      }
    }

    // ********* PRINTS **************:
    //cout << endl << "Counter = " << counter << "\t VS \t nb of pairs = " << k*(k-1)/2 << endl << endl;
    cout << "Best Merging: LogE = " << Diff_LogE_best << endl;
    cout << "Nb of pairs = " << k*(k-1)/2 << endl;
    /*
    for(it = Buffer_Community.begin(); it != Buffer_Community.end(); it++)
    {
      cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
    } */

    // ********* PERFORM MERGING *****:
    it = Buffer_Community.begin();  i_keep = (*it).first;   //  new_part = (*it).second;
    it++;  i_erase = (*it).first; //  new_part += (*it).second;

    Best_MCM_Partition[i_keep] = Best_MCM_Partition[i_keep] + Best_MCM_Partition[i_erase];
    Best_MCM_Partition.erase(i_erase);

    Best_MCM_LogE_Part[i_keep] = LogE_merged_best;
    Best_MCM_LogE_Part.erase(i_erase);

    // ********* STOPPING CRITERIA **************:
    if (Diff_LogE_best == 0) {  stop = true; }
    iteration++;
  } // End While

  // ********* FINAL PARTITION **************:
  cout << endl << "******** FINAL PARTITION -- iteration " << iteration << ":" << endl << endl;

  LogE_tot = 0;
  for(it = Best_MCM_Partition.begin(); it != Best_MCM_Partition.end(); it++)
  {
    cout << (*it).first << "\t " << bitset<n>((*it).second) << endl;
    LogE_tot += Best_MCM_LogE_Part[(*it).first];
  } 

  cout << endl << "\t Number of parts = " << Best_MCM_Partition.size() << endl;
  cout << "\t LogE = " << LogE_tot << endl;   // LogE_MCM(Kset, Best_MCM_Partition, N) 

}

/******************************************************************************/
/*******************************   main function   ****************************/
/******************************************************************************/
int main()
{  
  cout << "--->> Create OUTPUT Folder: (if needed) ";
  system( ("mkdir -p " + OUTPUT_directory).c_str() );
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
  unsigned int N=0; // will contain the number of datapoints in the dataset
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
  for (unsigned int i = 0 ; i<r; i++)
  {
    Partition_Indep[i] = Op;
    //cout << "Added Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;    
    Op = Op << 1;
  }
  cout << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "****************************************  Merging 2:  *************************************";
  cout << endl << "*******************************************************************************************" << endl;
  /*
  // ************ merging2:
  auto start = chrono::system_clock::now();
  Hierarchical_merging2(Kset, N);
  
  auto end = chrono::system_clock::now();
  chrono::duration<double> elapsed = end - start;
  cout << endl << "Elapsed time: " << elapsed.count() << "s";
*/
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "****************************************  Merging 3:  *************************************";
  cout << endl << "*******************************************************************************************" << endl;
  // ************ merging3:
  auto start = chrono::system_clock::now();
  Hierarchical_merging3(Kset, N);
  
  auto end = chrono::system_clock::now();
  chrono::duration<double> elapsed = end - start;
  cout << endl << "Elapsed time: " << elapsed.count() << "s";

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "****************************************  Merging 4:  *************************************";
  cout << endl << "*******************************************************************************************" << endl;
  // ************ merging3:
 /* auto start = chrono::system_clock::now();
  Hierarchical_merging4(Kset, N);
  
  auto end = chrono::system_clock::now();
  chrono::duration<double> elapsed = end - start;
  cout << endl << "Elapsed time: " << elapsed.count() << "s";
*/
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "********************************  Big 5, random communities:  *****************************";
  cout << endl << "*******************************************************************************************" << endl;

  initialise_generator();
//  unsigned int A = 5;

  map<unsigned int, uint64_t> Partition_rand, Partition_rand2;

/*
  Partition_rand = Random_Partition(A);
  for(map<unsigned int, uint64_t>::iterator it = Partition_rand.begin(); it != Partition_rand.end(); it++)
  {
    cout << (*it).second << "\t " << bitset<n>((*it).second) << endl;
  }
  cout << "Is this partition correct?  " << ((check_partition(Partition_rand))?"Yes":"No") << endl;
  cout << " \t LogE = " << LogE_MCM(Kset, Partition_rand, N) << " \t LogL = " << LogL_MCM(Kset, Partition_rand, N) << endl;

  Partition_rand.clear();
*/
/*
  double LogE_now=0, min_LogE = 0, max_LogE = 0; 

  fstream file("./OUTPUT/LogE_RandomPart.dat", ios::out);

  Partition_rand = Random_Partition();
  min_LogE = LogE_MCM(Kset, Partition_rand, N);
  max_LogE = min_LogE;

  file << min_LogE << endl;

  int x = 100;
  for (int i=0; i<(x-1); i++)
  {
    Partition_rand = Random_Partition();
    LogE_now = LogE_MCM(Kset, Partition_rand, N);   
    file << LogE_now << endl;

    if(LogE_now > max_LogE) { max_LogE = LogE_now; }
    if(LogE_now < min_LogE) { min_LogE = LogE_now; }

    Partition_rand.clear();
  }
  file.close();

  cout << "Number of sampled partitions = " << x << endl;
  cout << "Min(LogE) = " << min_LogE << endl;
  cout << "Max(LogE) = " << max_LogE << endl;
*/

/*
  Partition_rand2 = Random_Partition(A);
  for(map<unsigned int, uint64_t>::iterator it = Partition_rand2.begin(); it != Partition_rand2.end(); it++)
  {
    cout << (*it).second << "\t " << bitset<n>((*it).second) << endl;
  }
  cout << " \t LogE = " << LogE_MCM(Kset, Partition_rand2, N) << " \t LogL = " << LogL_MCM(Kset, Partition_rand2, N) << endl;
*/

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "**************************  All Successive Sub-Complete Models:  **************************";
  cout << endl << "*******************************************************************************************" << endl;
/*
  cout << "Sub-Complete models in the new basis:" << endl;

  map<uint32_t, uint32_t> Partition_SC;  Op = 1;
  Partition_SC[0] = 0;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_SC[0] += Op;
    cout << "Added Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_SC, N) << " \t LogL = " << LogL_MCM(Kset, Partition_SC, N) << endl;    
    Op = Op << 1;
  }
*/
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "********************************  Compare all MCModels:  **********************************";
  cout << endl << "*******************************************************************************************" << endl;
/*
// *** Search among all MCM based on the r first basis operators (models of rank exactly equal to `r`):
  int r=9;
  double LogE_BestMCM = 0;
  map<uint32_t, uint32_t> MCM_Partition = MCM_GivenRank_r(Kset, r, N, &LogE_BestMCM);
//cout << "\t Best LogE = " << LogE_BestMCM << endl;

// *** Search among all MCM based on the r first basis operators,
// ***  where   1 <= r <= m, where m is the size of Basis_li:
  //map<uint32_t, uint32_t> MCM_Partition = MCM_allRank(Kset, N, &LogE_BestMCM);

// *** 
  PrintTerminal_MCM_Info(Kset, N, MCM_Partition);
*/
  return 0;
}
