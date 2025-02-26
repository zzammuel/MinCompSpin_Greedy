#define _USE_MATH_DEFINES
#include <bitset>
#include <cmath>       /* tgamma */
#include <map>

using namespace std;

#include "data.h"


/******************************************************************************/
/**************************   MODEL COMPLEXITY   ******************************/
/******************************************************************************/
double GeomComplexity_SubCM(unsigned int m);

/******************************************************************************/
/**************** Log-Evidence (LogE) of a sub-complete model  ****************/
/******************************************************************************/
// Compute the log-evidence of a sub-complete model based on m basis elements
// ! Kset must have been previously reduced to these m basis elements !
// This function is mainly used for call by `LogE_PartMCM`,
// but can also be used to compute the log-likelihood of a complete model
//
double LogE_SubC_forMCM(map<__int128_t, unsigned int> Kset, uint32_t m, unsigned int N)
{
  double LogE = 0;

  map<__int128_t, unsigned int>::iterator it;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;

  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) {cout << "problem Ks = 0 for some mu_m" << endl; }
    LogE += lgamma(Ks + 0.5);
  }
  if (Ncontrol != N) { cout << "Error Likelihood function: Ncontrol != N" << endl;  }

//  return LogE - GeomComplexity_SubCM(m) - lgamma( (double)( N + (1UL << (m-1)) ) );
  return LogE + lgamma((double)( 1UL << (m-1) )) - (Kset.size()/2.) * log(M_PI) - lgamma( (double)( N + (1UL << (m-1)) ) ); 
}

/******************************************************************************/
/*********  Log-Evidence (LogE) of a sub-complete part of a MCM   *************/
/******************************************************************************/
// Compute the log-evidence of the sub-complete part (of an MCM) defined by Ai.
// This function could be also used directly by the user
// to compute the log-likelihood of a sub-complete model
// Rem: the function compute the LogE as if the space were reduced to the sub-space defined by the model

double LogE_SubCM(map<__int128_t, unsigned int> Kset, __int128_t Ai, unsigned int N, bool print_bool = false)
{
  map<__int128_t, unsigned int>::iterator it;
  map<__int128_t, unsigned int> Kset_new;

  __int128_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset

  if (print_bool)  { 
  cout << endl << "--->> Build Kset for SC Model based on " << bitset<n>(Ai) << " for MCM.." << endl;
  }
//Build Kset:
  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    s = it->first;      // initial state s 
    ks = it->second;    // # of times s appears in the data set
    if (print_bool)  {  cout << bitset<n>(s) << " \t" ;  }

    s &= Ai;   // troncated state: take only the bits indicated by Ai
//    sig_m = bitset<m>(bitset<m>(mu).to_string()).to_ulong(); //bitset<m>(mu).to_ulong(); // mu|m
    if (print_bool)  {  cout << bitset<n>(s) << endl; }

    Kset_new[s] += ks;
    //Kset[mu_m].second.push_back(make_pair(mu, N_mu));
  }
  if (print_bool)  {  cout << endl;  }

  bitset<n> hi{ static_cast<unsigned long long>(Ai >> 64) },
      lo{ static_cast<unsigned long long>(Ai) },
      bits{ (hi << 64) | lo };

  return LogE_SubC_forMCM(Kset_new, bits.count(), N);
}

/******************************************************************************/
/****************************   LogE of a MCM   *******************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the r elements, 
// i.e., that no basis element appears in more than 1 part of the partition.

unsigned int count_bits(__int128_t bool_nb)
{
  bitset<n> hi{ static_cast<unsigned long long>(bool_nb >> 64) },
            lo{ static_cast<unsigned long long>(bool_nb) },
            bits{ (hi << 64) | lo };
  return bits.count();
}
/*
pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition)
{
  map<unsigned int, __int128_t>::iterator Part;
  __int128_t sum = 0;
  unsigned int rank = 0; 

  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += count_bits((*Part).second); 
  }

  return make_pair((count_bits(sum) == rank), rank); 
}*/

pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition)
{
  map<unsigned int, __int128_t>::iterator Part;
  __int128_t sum = 0;
  unsigned int rank = 0; 

  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += count_bits((*Part).second); 
  }

  return make_pair((count_bits(sum) == rank), rank); 
}

double LogE_MCM(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r = n, bool print_bool = false)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogE = 0; 
    unsigned int rank = 0;
    map<unsigned int, __int128_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
        bitset<n> hi{ static_cast<unsigned long long>((*Part).second >> 64) },
            lo{ static_cast<unsigned long long>((*Part).second) },
            bits{ (hi << 64) | lo };

      LogE += LogE_SubCM(Kset, (*Part).second, N);
      rank += bits.count();
    }  
    return LogE - ((double) (N * (n-rank))) * log(2.);
  //}
  return 0;
}

