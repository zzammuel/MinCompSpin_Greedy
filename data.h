#include <iostream>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/

// number of binary (spin) variables:
const unsigned int n = 9; //120;  

// INPUT DATA FILES (optional):  
// the input datafile can also be specified directly in the main() function, as an argument of the function "read_datafile()":
const string datafilename = "INPUT/SCOTUS_n9_N895_Data.dat"; //"INPUT/sampled.dat";

// INPUT BASIS FILES (optional):
const string basis_IntegerRepresentation_filename = "INPUT/SCOTUS_n9_BestBasis_Integer.dat";        // (optional) Input basis file 
const string basis_BinaryRepresentation_filename = "INPUT/SCOTUS_n9_BestBasis_Binary.dat";      // (optional) Input basis file

// INPUT MCM FILES (optional): // For example: choice of a community to compare with the uncovered community:
const string communityfile = "INPUT/SCOTUS_Communities_inBestBasis.dat"; //"INPUT/community.dat";

// OUTPUT FOLDER:
const string OUTPUT_directory = "OUTPUT/";


/********************************************************************/
/**********************    DO NOT CHANGE    *************************/
/********************************************************************/
const __int128_t un = 1;
const __int128_t NOp_tot = (un << n) - 1;


/********************************************************************/
/*********************    Proba Structure    ************************/
/********************************************************************/
struct Proba {
    __int128_t s;         // state in the original basis
    __int128_t sig;       // state in the new basis

    double P_D_s = 1.;    // empirical probability of s  
    double P_D_sig = 1.;  // empirical probability of sig --> this should be the same then s if r=n, i.e. if the MCM models all the n spins
    double P_MCM = 1.;  // model probability of s 
};


/********************************************************************/
/****************    FOR METROPOLIS ALGORITHM    ********************/
/********************************************************************/
    struct Interaction
    {
        __int128_t Op;      // binary operator associated to the interaction
        double g;   // parameter of the interaction in {-1,+1} representation
        double av_M;      // average in the Model
        double av_D;      // average in the generated Data
    };


/********************************************************************/
/********************    NETWORK CONSTANTS    ***********************/
/********************************************************************/
// File to store test results
const string GNdatafile = "INPUT/Test.dat";

//const unsigned int nNodes = 20;

// Files for generating network data
const string networkfile = "INPUT/network.dat";


/********************************************************************/
/****************    FOR CHOICE OF THE BEST BASIS    ****************/
/********************************************************************/

    //Structure with the final information for the probability of appearance of each operator in the dataset
struct Operator
{
  uint32_t bin;     // binary representation of the operator
  mutable unsigned int layer;        // to which layer the operator belongs --> known only after the selection of the best Basis: by default, equal to n (=last layer)
  //unsigned int k1;  // nb of point where op = 1 --> it's a R.V.:  k1 = sum(op[s^i])
  mutable double p1_M;     // in the model: probability that op = 1 
  double p1_D;     // in the data: probability that op = 1 --> rem: it's a R.V. = sum(op[s^i]) / N

  double S;           // - [ p1*log(p1) + (1-p1)*log(1-p1) ]
  mutable double DKL;
  bool operator < (const Operator &other) const   // for ranking Operators from the most to the less likely
    { return S <= other.S; }
};

struct sort_by_prob
{
    bool operator()(const Operator& a, const Operator&b) const {
        return abs(1/2 - a.p1_D) >= abs(1/2 - b.p1_D);
    }
};

const unsigned int alpha = 3;

/********************************************************************/
/**********************    CURRENTLY NOT IN USE    ******************/
/********************************************************************/

    const string basis_FromIndices_filename = "INPUT/Big5PT_Best_Basis.dat";
    const string MCM_FromIndices_filename = "INPUT/Big5PT_Best_MCM_7.dat";

