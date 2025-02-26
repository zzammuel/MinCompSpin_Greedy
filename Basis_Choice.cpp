#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <bitset>

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
list<__int128_t> Read_BasisOp_BinaryRepresentation(string Basis_binary_filename = basis_BinaryRepresentation_filename)
{
  __int128_t Op_bit = 0, Op = 0;
  list<__int128_t> Basis_li;

  ifstream myfile (Basis_binary_filename.c_str());
  string line, line2;   
  char c = '1';  

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);          //take the n first characters of line

      //convert string line2 into a binary integer:
      Op_bit = un << (n - 1);
      Op = 0;
      for (auto &elem: line2)
      {
        if (elem == c) { Op += Op_bit; }
        Op_bit = Op_bit >> 1;
      }

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
list<__int128_t> Read_BasisOp_IntegerRepresentation(string Basis_integer_filename = basis_IntegerRepresentation_filename)
{
  __int128_t Op = 0;
  list<__int128_t> Basis_li;

  ifstream myfile (Basis_integer_filename.c_str());
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
/*************************    Original Basis     ******************************/
/******************************************************************************/
list<__int128_t> Original_Basis()
{
  __int128_t Op = 1;
  list<__int128_t> Basis_li;

  for (int i=0; i<n; i++)
  {
    Basis_li.push_back(Op);
    Op = Op << 1;
  }

  return Basis_li;
}

/******************************************************************************/
/***************************    Print Basis     *******************************/
/******************************************************************************/
void PrintTerm_Basis(list<__int128_t> Basis_li)
{
  int i = 1;
  for (list<__int128_t>::iterator it = Basis_li.begin(); it != Basis_li.end(); it++)
  {
    cout << "##\t " << i << " \t " << bitset<n>(*it) << endl; i++;
  } 
  cout << "##" << endl;
}

