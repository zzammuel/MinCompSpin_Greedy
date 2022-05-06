#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <list>
#include <bitset>
#include <map>
#include <cstring>

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

list<Interaction> write_interactions(double J, string file)
{
    string line, line2;
    list<Interaction> I_list;
    Interaction I;

    __int128_t Op = 1;
    Op <<= (n - 1);

    bool flag = false;

    // Open file to read network
    ifstream myfile(file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            I.Op = 0;
            I.g = J;
            while (getline(ss, line2, '\t'))
            {
                if ((Op >> stoi(line2) - 1) < I.Op) { flag = true; }
                I.Op += (Op >> (stoi(line2) - 1));
            }
            if (flag) { flag = false; continue; }
            I_list.push_back(I);
        }
        myfile.close();
    }
    return I_list;
}

map<unsigned int, __int128_t> read_communities(string file)
{
    map<unsigned int, __int128_t> Partition;

    string line, line2;
    __int128_t Op = 1;
    Op <<= n - 1;
    vector<int> comm;

    ifstream myfile(file.c_str());
    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            stringstream ss(line);
            while (getline(ss, line2, '\t'))
            {
                comm.push_back(stoi(line2));
            }
            Partition[comm[1]] += Op;
            Op >>= 1;

            comm.clear();
        }
        myfile.close();
    }
    return Partition;
}

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
map<__int128_t, unsigned int> read_datafile(unsigned int *N)    // O(N)  where N = data set size
{
    string line, line2;     char c = '1';
  __int128_t nb = 0, Op;
  (*N) = 0;            // N = dataset size
  //cout << endl << "--->> Read \"" << datafilename << "\",\t Build Nset...";

// ***** data are store in Nset:  ********************************
  map<__int128_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

  ifstream myfile (datafilename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
        Op = 1;
        Op = (Op << n - 1);
        nb = 0;
      line2 = line.substr (0,n);          //take the n first characters of line
      for (auto &elem: line2)
      {
          if (elem == c) { nb += Op; }
          Op = Op >> 1;
      }
      Nset[nb] += 1;
      //cout << line << endl;   //cout << nb << " :  " << bitset<n>(nb) << endl;
      (*N)++;
    }
    myfile.close();
  }
  else
  {
      cout << "Unable to open file";
      (*N) = 0;
  }

  //cout << "\t\t data size N = " << (*N) << endl;

  return Nset;
}

/*
void Print_File_Nset(map<__int128_t, unsigned int> Nset)
{
  map<__int128_t, unsigned int>::iterator it;
  it=Nset.begin();
  __int128_t s;

//  s=it->first;
//  cout << s << endl;  

//  for(it=Nset.begin(); it!=Nset.end(); it++)
  for (int i=0; i<10; i++)
  {
    //s=it->first;
    cout << (*it).first << "\t" << (*it).second << endl;
    it++;
  }
}*/

void Print_File_Nset (map<__int128_t, unsigned int> Nset, unsigned int N, string OUTPUTfilename)
// map.second = nb of time that the state map.first appears in the data set
{
  map<__int128_t, unsigned int>::iterator it;
  int Ncontrol = 0;
  __int128_t un = 1;

  fstream file(OUTPUTfilename.c_str(), ios::out);
  file << "#N = " << N << endl;
  file << "#Total number of accessible states = 2^(" << n << ") - 1" << endl;
  file << "#Number of visited states, Nset.size() = " << Nset.size() << endl;
  file << "#" << endl;
  file << "#1: state \t #2: nb of pts in state \t #3: Pba state" << endl;

  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    file << bitset<n>((*it).first) << " => " << (*it).second; // << endl;
    file << "  \t  P = " << ((*it).second) / (float) N << endl;
    Ncontrol += (*it).second;
  }

  if (Ncontrol != N) { cout << "Error function \'read_Nset\': Ncontrol != N" << endl;  }

  file.close();
}

/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a model (defined by the m basis vector) --> return the new m-state (state in the new m-basis)
// Rem: must have m <= n 
__int128_t transform_mu_basis(__int128_t mu, list<__int128_t> basis)
{
  __int128_t un_i = 1;
  __int128_t final_mu = 0;

  list<__int128_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    if ( (bitset<n>( (*phi_i) & mu ).count() % 2) == 1) // odd number of 1, i.e. sig_i = 1
      {
        final_mu += un_i;
      }
    un_i = (un_i << 1);
  }

  return final_mu;
}

/******************************************************************************/
/************************** K_SET *********************************************/
/******************************************************************************/
// Build Kset for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:

map<__int128_t, unsigned int> build_Kset(map<__int128_t, unsigned int> Nset, list<__int128_t> Basis, bool print_bool=false)
// sig_m = sig in the new basis and cut on the m first spins 
// Kset[sig_m] = #of time state mu_m appears in the data set
{
  map<__int128_t, unsigned int>::iterator it;
  map<__int128_t, unsigned int > Kset;

  __int128_t s;        // initial state
  __int128_t sig_m;    // transformed state and to the m first spins

  unsigned int ks=0; // number of time state s appear in the dataset

  cout << endl << "--->> Build Kset..." << endl;

//Build Kset:
  for (it = Nset.begin(); it!=Nset.end(); ++it)
  {
    s = it->first;       // state s
    ks = it->second;    // # of times s appears in the data set
    sig_m = transform_mu_basis(s, Basis);
//    sig_m = bitset<m>(bitset<m>(mu).to_string()).to_ulong(); //bitset<m>(mu).to_ulong(); // mu|m
    if (print_bool)  {  cout << bitset<n>(s) << " \t" << ": \t" << bitset<n>(sig_m) << endl; }

    Kset[sig_m] += ks;
    //Kset[mu_m].second.push_back(make_pair(mu, N_mu));
  }
  cout << endl;

  return Kset;
}
