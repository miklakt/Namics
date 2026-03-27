
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
#include <array>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <span>
#include <nlohmann/json.hpp>

#include <iomanip>
#include <memory>
#include "debug_log.h"

// #define LongReal

using ParameterStore = nlohmann::ordered_json;
//these are our options.
#ifdef LongReal
typedef long double Real; //See comment at top of file.
#else
typedef double Real;
#endif



//I.V. Ionova, E.A. Carter, "Error std::vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996.

extern Real e;
extern Real T;
extern Real k_B;
extern Real k_BT;
extern Real PIE;
extern Real eps0;
enum MoleculeType {monomer, linear, branched};
enum transfer {to_segment, to_bm, reset};
enum EngineType {SCF};
enum LatticeType {simple_cubic, hexagonal};
enum CP{co_solvent_theta,co_solvent_phibulk,chi_C_D};
#endif
