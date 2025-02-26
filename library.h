#include <set>



/******************************************************************************/
/******************************************************************************/
/*********************     SPECIFY the BASIS    *******************************/
/******************************************************************************/
/******************************************************************************/

/*** Original Basis:    ***********************************************/
/******************************************************************************/
list<__int128_t> Original_Basis();   // return the original basis, i.e., {s1, s2, ..., sn}

/*** READ BASIS from a FILE:    ***********************************************/
/******************************************************************************/
list<__int128_t> Read_BasisOp_BinaryRepresentation(string Basis_binary_filename = basis_BinaryRepresentation_filename);   // filename to specify in data.h
list<__int128_t> Read_BasisOp_IntegerRepresentation(string Basis_integer_filename = basis_IntegerRepresentation_filename); 

/*** Print Basis Info in the Terminal:    *************************************/
/******************************************************************************/
void PrintTerm_Basis(list<__int128_t> Basis_li);




/******************************************************************************/
/******************************************************************************/
/******************     READ and TRANSFORM DATA    ****************************/
/******************************************************************************/
/******************************************************************************/

/*** READ DATA and STORE data in Nset:    *************************************/
/******************************************************************************/
map<__int128_t, unsigned int> read_datafile(unsigned int *N, string file = datafilename); // filename to specify in data.h

/*** DATA CHANGE of BASIS:    *************************************************/
/******************************************************************************/
// *** Build Kset with the following definitions:
// *** mu_m = states of the systems written in the basis specified in `list<uint32_t> Basis`
// *** Kset[sig_m] = Number of times the state mu_m appears in the transformed dataset
//
// *** Rem: the new basis can have a lower dimension then the original dataset; 
// *** in which case the function will reduce the dataset to the subspace defined by the specified basis.
map<__int128_t, unsigned int> build_Kset(map<__int128_t, unsigned int> Nset, list<__int128_t> Basis, bool print_bool=false);




/******************************************************************************/
/******************************************************************************/
/***************** Log-LIKELIHOOD (LogL), Log-EVIDENCE (LogE) *****************/
/***************************  and COMPLEXITY   ********************************/
/******************************************************************************/
/******************************************************************************/

/****************   for a sub-Complete Model (SubCM)   ************************/
/**********  restricted to the subspace of the Sub-Complete Model  ************/
/******************************************************************************/
// *** the SubCM is the one specified in Ai;
// *** Ai must be an integer encoded on at least n bits, where each 1 indicates the basis elements included in the part:
// *** For ex. Ai = 01001 is encoded on n=5 basis elements, and element Op1 and Op4 belong to the part;
// *** Rem: Basis elements are ordered from the right to the left.

double LogL_SubCM(map<__int128_t, unsigned int > Kset, __int128_t Ai, unsigned int N, bool print_bool = false);
double LogE_SubCM(map<__int128_t, unsigned int > Kset, __int128_t Ai, unsigned int N, bool print_bool = false);

// *** Complexity of a SC model based on m basis Operators: m >= 1. Rem: C_geom(m=1) = log(pi):
double GeomComplexity_SubCM(unsigned int m);                  // Geometric complexity
double ParamComplexity_SubCM(unsigned int m, unsigned int N); // Complexity due to the number of parameters

/******************   for a Complete Model (CM)   *****************************/
/******************************************************************************/
double LogL_CM(map<__int128_t, unsigned int > Kset, unsigned int N);

/****************************    for a MCM     ********************************/
/******************************************************************************/
double LogL_MCM(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r = n, bool print_bool = false);
double LogE_MCM(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, unsigned int r = n, bool print_bool = false);
double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, double *C_param, double *C_geom);




/******************************************************************************/
/******************************************************************************/
/********************   DEFINE MCMs and PRINT INFO   **************************/
/******************************************************************************/
/******************************************************************************/

unsigned int count_bits(__int128_t bool_nb);

// *** Define an MCM from a file; Each part must be encoded in a binary number over n spins:
//map<uint32_t, uint32_t> Read_MCMParts_BinaryRepresentation(string MCM_binary_filename);

// *** Check that the provided model corresponds to a partition of the basis variables (i.e. properly defines an MCM):
pair<bool, unsigned int> check_partition(map<unsigned int, __int128_t> Partition);  // the second element is the rank of the partition (dimension of the MCM)

// *** Print information about the MCM specified in `MCM_Partition`:
void PrintTerminal_MCM_Info(map<__int128_t, unsigned int> Kset, unsigned int N, map<unsigned int, __int128_t> MCM_Partition);
void Print_Partition(map<unsigned int, __int128_t> partition);

// *** Create successive independent models defined on the new basis, and print the corresponding information:
//void PrintInfo_All_Indep_Models(map<uint32_t, unsigned int> Kset, unsigned int N);

// *** Create successive Sub-complete models defined on the new basis, and print the corresponding information:
//void PrintInfo_All_SubComplete_Models(map<uint32_t, unsigned int> Kset, unsigned int N);

map<unsigned int, __int128_t> read_communities(string file);

// Check if subset
bool is_subset(map<unsigned int, __int128_t> fp1, map<unsigned int, __int128_t> fp2);


/******************************************************************************/
/******************************************************************************/
/***************************   Find Best MCM   ********************************/
/******************************************************************************/
/******************************************************************************/

map<unsigned int, __int128_t> Matrix(map<__int128_t, unsigned int> Kset, unsigned int N, unsigned int r = n);



/******************************************************************************/
/******************************************************************************/
/********************   Statistical Properties   ******************************/
/******************************************************************************/
/******************************************************************************/

// *** Statistical properties
double Entropy(map <__int128_t, unsigned int> Kset, unsigned int N);
double Kullback_Leibler(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N);
double JS_divergence(map<__int128_t, double> Prob1, map<__int128_t, double> Prob2, unsigned int N);

double Var_of_Inf(map<unsigned int, __int128_t> Partition1, map<unsigned int, __int128_t> Partition2);
double Norm_Mut_info(map<unsigned int, __int128_t> Partition1, map<unsigned int, __int128_t> Partition2);



/******************************************************************************/
/******************************************************************************/
/***************************   PRINT TO FILE:  ********************************/
/******************   DATA VS MODEL STATE PROBABILITIES  **********************/
/******************************************************************************/
/******************************************************************************/
void PrintFile_MCM_Info(list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, string filename = "Result");

void PrintFile_StateProbabilites_NewBasis(map<__int128_t, unsigned int > Kset, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, string filename = "Result");
void PrintFile_StateProbabilites_OriginalBasis(map<__int128_t, unsigned int > Nset, list<__int128_t> Basis, map<unsigned int, __int128_t> MCM_Partition, unsigned int N, string filename = "Result");

// *** Create distributions of empirical data and MCMs
map<__int128_t, double> emp_dist(map<__int128_t, unsigned int> Kset, unsigned int N);
map<__int128_t, double> MCM_distr(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N);




/******************************************************************************/
/******************************************************************************/
/**********************   Metropolis algorithm   ******************************/
/******************************************************************************/
/******************************************************************************/

list<Interaction> write_interactions(double J, string file);

