
#ifndef NAMICSxH
#define NAMICSxH
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <span>

#include <iomanip>
#include "debug_log.h"

using namespace std;
using ParameterStore = std::unordered_map<std::string, std::string>;
//these are our options.
#ifdef LongReal
typedef long double Real; //See comment at top of file.
#else
typedef double Real;
#endif



//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996.

#ifndef MAINxH //here define global variables.
extern Real* BlasResult;
extern string version;
extern Real e;
extern Real T;
extern Real k_B;
extern Real k_BT;
extern Real PIE;
extern int DEBUG_BREAK;
extern Real eps0;
extern bool debug;
#endif

enum MoleculeType {monomer, linear, branched, dendrimer, asym_dendrimer, comb, water};
enum transfer {to_segment, to_bm, reset};
enum EngineType {SCF, MICRO};
enum LatticeType {simple_cubic, hexagonal};
enum SolverType {Pseudohessian,Lbfgs};
enum CP{co_solvent_theta,co_solvent_phibulk,chi_C_D};
#endif
