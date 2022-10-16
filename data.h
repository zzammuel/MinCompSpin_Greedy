#include <iostream>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
const unsigned int n = 120;                          // number of spin variables
//const unsigned int nNodes = 20;

const __int128_t un = 1;
const __int128_t NOp_tot = (un << n) - 1;


// File to store test results
const string GNdatafile = "INPUT/Test.dat";

    // Input datafile
    const string datafilename = "INPUT/sampled.dat";

    const string OUTPUT_directory = "INPUT/";

    // Files for generating network data
    const string networkfile = "INPUT/network.dat";

    // Exact community
    const string communityfile = "INPUT/community.dat";

    struct Interaction
    {
        __int128_t Op;      // binary operator associated to the interaction
        double g;   // parameter of the interaction in {-1,+1} representation
        double av_M;      // average in the Model
        double av_D;      // average in the generated Data
    };

    // Not used
    const string basis_IntegerRepresentation_filename = "INPUT/SCOTUS_n9_Basis_Integer.dat";		// Input basis file
    const string basis_BinaryRepresentation_filename = "INPUT/SCOTUS_n9_Basis_Binary.dat";		// Input basis file

    const string basis_FromIndices_filename = "INPUT/Big5PT_Best_Basis.dat";
    const string MCM_FromIndices_filename = "INPUT/Big5PT_Best_MCM_7.dat";

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