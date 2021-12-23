#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <bitset>
#include <map>

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/*****************   Read Basis Operators from file  **************************/
/******************************************************************************/

/******  1)  Operators should be written in one single column          ********/
/******  2)  Operators can be written in two different versions:       ********/
/***   a) as a binary representation of the spin involved in the Operator;    */
/***   b) as the integer value of that binary representation.                 */
/******************************************************************************/
/****  Ex. for a system with n=4 spin:  ***************************************/
/****      -- the operator s1 would be encoded as 0001,                      **/
/****      which has the integer representation 1  -->   0001 = 1      ********/
/****      -- the operator s1 s2 would be encoded as 0011,             ********/
/****      which has the integer representation 3  -->   0011 = 3      ********/
/******************************************************************************/


/******************************************************************************/
/*** VERSION a) Operators are written as the binary          ******************/
/****           representation of the interactions           ******************/
/******************************************************************************/
list<uint32_t> Read_BasisOp_BinaryRepresentation()
{
  uint32_t Op = 0;
  list<uint32_t> Basis_li;

  ifstream myfile (basis_BinaryRepresentation_filename.c_str());
  string line, line2;     

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);          //take the n first characters of line

      Op = bitset<n>(line2).to_ulong();   //convert string line2 into a binary integer
      Basis_li.push_back(Op);   
    }
    myfile.close();
  }

  return Basis_li;
}

/******************************************************************************/
/*** VERSION b) Operators are written as the integer values of the binary *****/
/****           representation of the interactions           ******************/
/******************************************************************************/
list<uint32_t> Read_BasisOp_IntegerRepresentation()
{
  uint32_t Op = 0;
  list<uint32_t> Basis_li;

  ifstream myfile (basis_IntegerRepresentation_filename.c_str());
  string line;    

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      Op = stoi(line);
      Basis_li.push_back(Op);
    }
    myfile.close();
  }

  return Basis_li;
}

/******************************************************************************/
/*** VERSION c) Operators are written as the integer values of the binary *****/
/****           representation of the interactions           ******************/
/******************************************************************************/
list<uint64_t> Read_BasisOp_FromIndices()
{
  string line;
  int temp;    
  stringstream sstream; 

  uint64_t un = 1;
  uint64_t Basis_Op = 0;
  list<uint64_t> BasisVect_li;

  ifstream myfile (basis_FromIndices_filename.c_str()); 

  //cout << "Init: " << bitset<n>(un) << "  " << bitset<n>((un << (n-1))) << endl;
  //cout << "Init: " << bitset<n>(un) << "  " << bitset<n>((un << 0)) << endl;

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      sstream = stringstream(line);  // Storing the whole string into string stream
      Basis_Op = 0;
      while (sstream >> temp)   // extracting first "word" from stream && Checking the given "word" is integer or not
      { 
        Basis_Op += (un << (n-temp));
      }
      if(Basis_Op != 0)
      {
        BasisVect_li.push_back(Basis_Op);
        cout << Basis_Op << "\t" << bitset<n>(Basis_Op) << endl;
        //temp = "";
        //sstream.str("");
      }
    }
  }
  myfile.close();
  cout << endl;

  return BasisVect_li;
}

/******************************************************************************/
/*************************    Original Basis     ******************************/
/******************************************************************************/
list<uint64_t> Original_Basis()
{
  uint64_t Op = 1;
  list<uint64_t> Basis_li;

  for (int i=0; i<n; i++)
  {
    Basis_li.push_back(Op);
    Op = Op << 1;
  }

  return Basis_li;
}

/******************************************************************************/
/****************************    Print Basis     ******************************/
/******************************************************************************/
void PrintTerm_Basis(list<uint64_t> Basis_li)
{
  int i = 1;
  for (list<uint64_t>::iterator it = Basis_li.begin(); it != Basis_li.end(); it++)
  {
    cout << "##\t " << i << " \t " << (*it) << " \t " << bitset<n>(*it) << endl; i++;
  } cout << "##" << endl;
}

