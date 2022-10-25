#include <bitset>
#include <cmath>       /* tgamma */
#include <map>
#include <list>
#include <fstream>

using namespace std;

#include "data.h"


/******************************************************************************/
/****************   Return Kset over a chosen sub-basis b_a    ****************/
/******************************************************************************/
map<__int128_t, unsigned int> Build_Kset_ba(map<__int128_t, unsigned int> Kset, __int128_t Ai) // Ai = integer indicated the spin elements included in b_a
{
  map<__int128_t, unsigned int>::iterator it;
  map<__int128_t, unsigned int> Kset_new;

  __int128_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset

  cout << "--->> Build Kset for SC Model based on " << bitset<n>(Ai) << " for MCM.." << endl << endl;

//Build Kset_ba:
  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    s = it->first;      // initial state s 
    ks = it->second;    // # of times s appears in the data set

    s &= Ai;   // troncated state: take only the bits indicated by Ai
    Kset_new[s] += ks;
  }

  return Kset_new;
}

/******************************************************************************/
/****************************      Check partition     ************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition);


/******************************************************************************/
/*****************   Compute the contribution to P_MCM(s)   *******************/
/******************  due to the sub-CM defined by Kset_ba   *******************/
/******************************************************************************/
void update_proba_MCM(map<__int128_t, Proba> &all_P, map<__int128_t, unsigned int> Kset_ba, __int128_t Ai, unsigned int N)
{
  map<__int128_t, Proba>::iterator it_P;

  __int128_t s, sig;        // states
  unsigned int ks=0;      // number of time state s appear in the dataset
  double Nd = (double) N;

//Contribution of the basis "ba" of spins to the model probability P_MCM(s):
  for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)
  {
    s = it_P->first;      // initial state s 
    sig = s & Ai;         // troncated state: take only the bits indicated by Ai
    all_P[s].P_MCM *= Kset_ba[sig]/Nd;
  }
}

/******************************************************************************/
/******************************************************************************/
/*************      Compute and Print the model probability     ***************/
/*********************  for the MCM constructed from Kset  ********************/
/************************   with the given partition   ************************/
/******************************************************************************/
/******************************************************************************/

/*************      Compute the model probabilities     ***************/

// This function can be used directly on the original basis, by replacing Kset by Nset:
map<__int128_t, Proba> P_sig(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N) // Probabilities in the sigma basis
{
  // Fill in the data probability:
  map<__int128_t, Proba> all_P;
  map<__int128_t, unsigned int>::iterator it;
  double Nd = (double) N;

  __int128_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset

  // Check partition:
  //bool Is_partition = check_partition(Partition);
  pair<bool, unsigned int> Is_partition = check_partition(Partition);
  unsigned int rank = Is_partition.second;

  if (!Is_partition.first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    double pre_factor = 1./((double) (un << (n-rank))); 

    for (it = Kset.begin(); it!=Kset.end(); ++it)
    {   
      s = it->first;      // initial state s 
      ks = it->second;    // # of times s appears in the data set

      all_P[s].P_D_s = ((double) ks)/Nd;
      all_P[s].P_MCM = pre_factor;
    }

    // Compute the Kset over each part: Kset_ba:
    map<__int128_t, unsigned int> Kset_ba;
    map<unsigned int, __int128_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      Kset_ba = Build_Kset_ba(Kset, (*Part).second);         // (*Part).second) = Ai = integer indicated the spin elements included in b_a
      update_proba_MCM(all_P, Kset_ba, (*Part).second, N);
      Kset_ba.clear();
    }  
  }

  return all_P;
}

/*************      Print the model probabilities     ***************/

void PrintFile_StateProbabilites_NewBasis(map<__int128_t, unsigned int > Kset, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, string filename = "Result")
{
  // Probabilities in the sigma basis:
  map<__int128_t, Proba> P_all = P_sig(Kset, MCM_Partition, N);
  map<__int128_t, Proba>::iterator it_P;

  string Psig_filename = filename + "_DataVSMCM_Psig.dat";

  fstream file_P_sig((OUTPUT_directory + Psig_filename), ios::out);
  file_P_sig << "## 1:sig \t 2:P_D(sig) \t 3:P_MCM(sig)" << endl;

  for (it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    file_P_sig << bitset<n>(it_P->first) << "\t " << (it_P->second).P_D_s << "\t " << (it_P->second).P_MCM << endl;
  }

  file_P_sig.close();
}

/******************************************************************************/
/******************************************************************************/
/*******************      Compute the model probability     *******************/
/**************************  for the MCM constructed   ************************/
/***************   in a given basis, with a given partition   *****************/
/******************************************************************************/
/******************************************************************************/

__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis);


map<__int128_t, Proba> P_s(map<__int128_t, unsigned int> Nset, list<__int128_t> Basis, map<unsigned int, __int128_t> Partition, unsigned int N) // Probabilities in the sigma basis
{
  double Nd = (double) N;

  // Fill in the data probability:
  map<__int128_t, Proba> all_P;

  if (!check_partition(Partition).first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    // Build Kset:
    map<__int128_t, unsigned int > Kset;
    __int128_t s;        // initial state
    __int128_t sig_m;    // transformed state and to the m first spins
    unsigned int ks=0; // number of time state s appear in the dataset

  //Build Kset and fill in P[s] from the data:
    cout << "--->> Build Kset and fill in P[s] from the data..." << endl;

    map<__int128_t, unsigned int>::iterator it;
    for (it = Nset.begin(); it!=Nset.end(); ++it)
    {
      s = it->first;       // state s
      ks = it->second;     // # of times s appears in the data set
      sig_m = transform_mu_basis(s, Basis);

      // Fill in Kset:
      Kset[sig_m] += ks;

      // Fill in P[s] empirical and value of transformed state:
      all_P[s].P_D_s = ((double) ks)/Nd;
      all_P[s].sig = sig_m;  
    }

  // Compute the model probability of the MCM based on Kset using "Partition":
    cout << "--->> Compute P[s] for the chosen MCM..." << endl << endl;
    map<__int128_t, Proba> all_P_sig = P_sig(Kset, Partition, N);   // Compute P[s] for the MCM in the new basis

  // Report the values of P_MCM in the original all_P:
    map<__int128_t, Proba>::iterator it_P;
    for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)           // Report P[s] for the MCM in the original basis
    {
      sig_m = (it_P->second).sig;
      (it_P->second).P_MCM = all_P_sig[sig_m].P_MCM;
      (it_P->second).P_D_sig = all_P_sig[sig_m].P_D_s;  // not necessary a probability distribution anymore
    }
  
    all_P_sig.clear();
    Kset.clear();
  }

  return all_P;  
}


/******************************************************************************/
/*****************      PRINT FILE: INFO about an MCM     *********************/
/******************************************************************************/

void PrintFile_MCM_Info(list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, string filename = "Result")
{
  //***** PRINT BASIS: 
  fstream file_MCM_info((OUTPUT_directory + filename + "_MCM_info.dat"), ios::out);

  file_MCM_info << "## sig_vec = states in the chosen new basis (ideally the best basis), defined by the basis operators:" << endl;
  int i = 1;
  for (list<__int128_t>::iterator it = Basis.begin(); it != Basis.end(); it++)
  {
    file_MCM_info << "##\t sig_" << i << " = " << bitset<n>(*it) << endl; i++;
  } file_MCM_info << "##" << endl;

  // Print info about the model -- Print MCM:
  file_MCM_info << "## The MCM Partition is defined on the sig_vec basis by the following Parts:" << endl;

  //***** PRINT MCM: 
  i = 1;
  for (map<unsigned int, __int128_t>::iterator it = MCM_Partition.begin(); it != MCM_Partition.end(); it++)
  {    
    __int128_t Part = (*it).second;
    file_MCM_info << "##\t MCM_Part_" << i << " = " << bitset<n>(Part) << endl; i++;
  }
  file_MCM_info << "##" << endl;

  file_MCM_info.close();
}


/******************************************************************************/
/*************      Print the model probabilities in a file     ***************/
/******************************************************************************/
unsigned int count_bits(__int128_t bool_nb);

void PrintFile_StateProbabilites_OriginalBasis(map<__int128_t, unsigned int > Nset, list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, string filename = "Result")
{
  // Compute all the state probabilities:
  map<__int128_t, Proba> P_all = P_s(Nset, Basis, MCM_Partition, N);

  double *Pk_D = (double *)malloc((n+1)*sizeof(double)); 
  double *Pk_MCM = (double *)malloc((n+1)*sizeof(double)); 

  unsigned int k = 0;
  for(k=0; k<=n; k++)
  {
    Pk_D[k] = 0;
    Pk_MCM[k] = 0;
  }

  string Ps_filename = filename + "_DataVSMCM_Ps.dat";
  string Pk_filename = filename + "_DataVSMCM_Pk.dat";

  cout << "--->> Print information about the MCM in the file: \'" << filename << "_MCM_info.dat\'" << endl;
  cout << "--->> Print the state probabilities P(s) in the file: \'" << Ps_filename << "\'" << endl;
  cout << "--->> Print the probability of a state with k \'+1\' bits: \'" << Pk_filename << "\'" << endl << endl;

  //***** Print info about the model -- Print Basis and MCM:  **************/
  PrintFile_MCM_Info(Basis, MCM_Partition, filename);

  //***** Print P(s):  *****************************************************/
  __int128_t s;
  fstream file_Ps((OUTPUT_directory + Ps_filename), ios::out);

  file_Ps << "## s = states in the original basis" << endl;
  file_Ps << "## sig = states in the chosen new basis (ideally the best basis)" << endl;
  file_Ps << "## The chosen \'sig\'-basis and the chosen MCM are printed in the file " << filename << "_MCM_info.dat" << endl;

  // Print P(s): 
  file_Ps << "## " << endl;
  file_Ps << "## 1:s \t 2:P_D(s) \t 3:P_MCM(s) \t 4:sig" << endl;

  for (map<__int128_t, Proba>::iterator it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    s = it_P->first;
    file_Ps << bitset<n>(s) << "\t" <<  (it_P->second).P_D_s << "\t" << (it_P->second).P_MCM << "\t" << bitset<n>((it_P->second).sig) << endl;

    k = count_bits(s); //bitset<n>(s).count();
    Pk_D[k] += (it_P->second).P_D_s;      // P[k] in the data
    Pk_MCM[k] += (it_P->second).P_MCM;    // P[k] from the MCM
  }
  file_Ps.close();

  //***** Print P(k):   ***************************************************/
  fstream file_Pk((OUTPUT_directory + Pk_filename), ios::out);

  file_Pk << "## 1:k \t 2:P_D(k) \t 3:P_MCM(k)" << endl;

  for(k=0; k<=n; k++)
  {
    file_Pk << k << "\t" << Pk_D[k] << "\t" << Pk_MCM[k] << endl;
  }
  file_Pk.close();
}
