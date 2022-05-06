#include <iostream>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
const unsigned int n = 92;                          // number of spin variables
//const unsigned int nNodes = 20;

const __int128_t un = 1;
const __int128_t NOp_tot = (un << n) - 1;


// File to store test results
const string GNdatafile = "INPUT/GNtest.dat";

    // Input datafile
    const string datafilename = "INPUT/brain.dat";

    const string OUTPUT_directory = "OUTPUT/";

    // Files for generating network data
    const string networkfile = "network.dat";

    // Exact community
    const string communityfile = "community.dat";

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