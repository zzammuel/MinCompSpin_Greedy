
/******************************************************************************/
/******************     READ and TRANSFORM DATA    ****************************/
/******************************************************************************/

/*** READ DATA and STORE data in Nset:    *************************************/
/******************************************************************************/
map<uint64_t, unsigned int> read_datafile(unsigned int *N); // filename to specify in data.h

/*** READ BASIS from a FILE:    ***********************************************/
/******************************************************************************/
list<uint32_t> Read_BasisOp_BinaryRepresentation();   // filename to specify in data.h
list<uint32_t> Read_BasisOp_IntegerRepresentation(); 

/*** Original Basis:    ***********************************************/
/******************************************************************************/
list<uint64_t> Original_Basis();   // return the original basis, i.e., {s1, s2, ..., sn}

/*** Print Basis Info in the Terminal:    *************************************/
/******************************************************************************/
void PrintTerm_Basis(list<uint64_t> Basis_li);

/*** DATA CHANGE of BASIS:    *************************************************/
/******************************************************************************/
// *** Build Kset with the following definitions:
// *** mu_m = states of the systems written in the basis specified in `list<uint32_t> Basis`
// *** Kset[sig_m] = Number of times the state mu_m appears in the transformed dataset
//
// *** Rem: the new basis can have a lower dimension then the original dataset; 
// *** in which case the function will reduce the dataset to the subspace defined by the specified basis.
map<uint64_t, unsigned int> build_Kset(map<uint64_t, unsigned int> Nset, list<uint64_t> Basis, bool print_bool=false);

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

double LogL_SubCM(map<uint64_t, unsigned int > Kset, uint64_t Ai, unsigned int N, bool print_bool = false);
double LogE_SubCM(map<uint64_t, unsigned int > Kset, uint64_t Ai, unsigned int N, bool print_bool = false);

// *** Complexity of a SC model based on m basis Operators: m >= 1. Rem: C_geom(m=1) = log(pi):
double GeomComplexity_SubCM(unsigned int m);                  // Geometric complexity
double ParamComplexity_SubCM(unsigned int m, unsigned int N); // Complexity due to the number of parameters

/******************   for a Complete Model (CM)   *****************************/
/******************************************************************************/
double LogL_CM(map<uint64_t, unsigned int > Kset, unsigned int N);

/****************************    for a MCM     ********************************/
/******************************************************************************/
double LogL_MCM(map<uint64_t, unsigned int> Kset, map<unsigned int, uint64_t> Partition, unsigned int N, bool print_bool = false);
double LogE_MCM(map<uint64_t, unsigned int> Kset, map<unsigned int, uint64_t> Partition, unsigned int N, bool print_bool = false);
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
void PrintTerminal_MCM_Info(map<uint64_t, unsigned int> Kset, unsigned int N, map<unsigned int, uint64_t> MCM_Partition);

bool check_partition(map<unsigned int, uint64_t> Partition);





