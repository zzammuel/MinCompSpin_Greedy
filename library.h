
/******************************************************************************/
/******************     READ and TRANSFORM DATA    ****************************/
/******************************************************************************/

/*** READ DATA and STORE data in Nset:    *************************************/
/******************************************************************************/
map<__int128_t, unsigned int> read_datafile(unsigned int *N); // filename to specify in data.h

/*** READ BASIS from a FILE:    ***********************************************/
/******************************************************************************/
list<uint32_t> Read_BasisOp_BinaryRepresentation();   // filename to specify in data.h
list<uint32_t> Read_BasisOp_IntegerRepresentation(); 

/*** Original Basis:    ***********************************************/
/******************************************************************************/
list<__int128_t> Original_Basis();   // return the original basis, i.e., {s1, s2, ..., sn}

/*** Print Basis Info in the Terminal:    *************************************/
/******************************************************************************/
void PrintTerm_Basis(list<__int128_t> Basis_li);

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
/***************** Log-LIKELIHOOD (LogL), Log-EVIDENCE (LogE) *****************/
/***************************  and COMPLEXITY   ********************************/
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
double LogL_MCM(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, bool print_bool = false);
double LogE_MCM(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N, bool print_bool = false);
double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, double *C_param, double *C_geom);

/******************************************************************************/
/******************************   Find Best MCM   *****************************/
/******************************************************************************/
// *** Compute all partitions of a set using Algorithm H:

// *** Version 1: Compare all the MCM of rank r, 
// ***            based on the r first elements of the basis used to build Kset:
//map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int r, unsigned int N, double *LogE_best);

// *** Version 2: Compare all the MCM 
// ***            based on the r first elements of the basis used to build Kset
// ***            for all r=1 to basis.size()  
//map<uint32_t, uint32_t> MCM_allRank(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best);

// *** Print information about the MCM specified in `MCM_Partition`:
void PrintTerminal_MCM_Info(map<__int128_t, unsigned int> Kset, unsigned int N, map<unsigned int, __int128_t> MCM_Partition);

bool check_partition(map<unsigned int, __int128_t> Partition);


// *** Statistical properties
double Entropy(map <__int128_t, unsigned int> Kset, unsigned int N);
double Kullback_Leibler(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N);
double JS_divergence(map<__int128_t, double> Prob1, map<__int128_t, double> Prob2, unsigned int N);

double Var_of_Inf(map<unsigned int, __int128_t> Partition1, map<unsigned int, __int128_t> Partition2);
double Norm_Mut_info(map<unsigned int, __int128_t> Partition1, map<unsigned int, __int128_t> Partition2);

// *** Create distributions of empirical data and MCMs
map<__int128_t, double> emp_dist(map<__int128_t, unsigned int> Kset, unsigned int N);
map<__int128_t, double> MCM_distr(map<__int128_t, unsigned int> Kset, map<unsigned int, __int128_t> Partition, unsigned int N);


map<unsigned int, __int128_t> read_communities(string file);
list<Interaction> write_interactions(double J, string file);

// Check if subset
bool is_subset(map<unsigned int, __int128_t> fp1, map<unsigned int, __int128_t> fp2);

/********************************************************************/
/*********************    METROPOLIS.CPP    *************************/
/********************************************************************/
/*******************    SAMPLES A DATASET    ************************/
/********************************************************************/
map<unsigned int, list<Interaction>> write_interactions_metropolis(double J, string file);
double delta_energy(__int128_t state, list<Interaction> edges);
void sample_data_metropolis(double J, string input_file, string output_filename, unsigned int N = 1000);

/********************************************************************/
/*****************    GENERATE_DATA_EXACT.CPP    ********************/
/********************************************************************/
/*******************    SAMPLES A DATASET    ************************/
/********************************************************************/
// *** list_I = list of operators and parameters of the model
// *** N = nb of datapoints

// *** Samples a dataset with "N" datapoints:  ******************** //
void Sample_dataset(list<Interaction> list_I, string output_filename, unsigned int N = 1000);

// **************************************************************** //
// *** Samples a dataset with "N" datapoints,  ******************** //
// *** and compute the model averages and data averages  ********** //
// *** of all the model operators:                       ********** //
// *** fill this information in "list_I"                 ********** //
// **************************************************************** //
void Sample_dataset_AND_Print_ModelData_Info(list<Interaction>& list_I, string output_filename, unsigned int N = 1000);

/********************************************************************/
/************    PRINT INFORMATION ABOUT THE MODEL    ***************/
/**********    and the GENERATED DATA if there were    **************/
/********************************************************************/
double* Probability_AllStates_Ising(list<Interaction> list_I, double* Z);
int Op_Ising(uint32_t Op, uint32_t state);
void Model_averages_Ising(list<Interaction>& list_I);
