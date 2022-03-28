#include <iostream>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
const unsigned int n = 20;                          // number of spin variables
//const unsigned int nNodes = 20;

const uint32_t un = 1;
const uint32_t NOp_tot = (un << n) - 1;

const string testsfilename = "INPUT/vary_temp.dat";
const string GNdatafile = "INPUT/GNtest.dat";

    const string datafilename = "INPUT/benchmark.dat";   // SCOTUS_n9_N895_Data.dat";	// Input datafile
    const string basis_IntegerRepresentation_filename = "INPUT/SCOTUS_n9_Basis_Integer.dat";		// Input basis file
    const string basis_BinaryRepresentation_filename = "INPUT/SCOTUS_n9_Basis_Binary.dat";		// Input basis file

    const string basis_FromIndices_filename = "INPUT/Big5PT_Best_Basis.dat";
    const string MCM_FromIndices_filename = "INPUT/Big5PT_Best_MCM_7.dat";

    const string OUTPUT_directory = "OUTPUT/";

    // Files for generating network data
    const string networkfile = "network.dat";

    // Exact community
    const string communityfile = "community.dat";

    struct Interaction
    {
        uint32_t Op;      // binary operator associated to the interaction
        double g;   // parameter of the interaction in {-1,+1} representation
        double av_M;      // average in the Model
        double av_D;      // average in the generated Data
    };