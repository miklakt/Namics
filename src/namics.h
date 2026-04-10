
#ifndef NAMICSxH
#define NAMICSxH
#include <cmath>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <iterator>
#include <algorithm>
#include <span>
#include <nlohmann/json.hpp>

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
enum LatticeType {simple_cubic, hexagonal};
#endif
