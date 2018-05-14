#include "mesodyn.h"

//Constructor
Mesodyn::Mesodyn(vector<Input*> In_, vector<Lattice*> Lat_, vector<Segment*> Seg_, vector<Molecule*> Mol_, vector<System*> Sys_, vector<Newton*> New_, string name_) {
  In = In_;
  name = name_;
  Lat = Lat_;
  Mol = Mol_;
  Seg = Seg_;
  Sys = Sys_;
  New = New_;
  KEYS.push_back("timesteps");
  KEYS.push_back("timebetweensaves");
  xNeighbors.reserve(2);
  yNeighbors.reserve(2);
  zNeighbors.reserve(2);
  componentNo = 0;
}
Mesodyn::~Mesodyn() {
}

void Mesodyn::AllocateMemory() {
  if (debug)
    cout << "nothing to allocate in Mesodyn" << endl;
}

bool Mesodyn::CheckInput(int start) {
  if (debug)
    cout << "Check Mesodyn" << endl;
  bool success = true;
  //right now all the checks for engine are done in input.cpp. could be move back here.
  success = In[0]->CheckParameters("mesodyn", name, start, KEYS, PARAMETERS, VALUES);
  if (success) {
    vector<string> options;
    if (GetValue("timesteps").size() > 0) {
      success = In[0]->Get_int(GetValue("timesteps"), timesteps, 1, 10000, "The number of timesteps should be between 1 and 10000");
    }
    cout << "timesteps is " << timesteps << endl;
  }
  return success;
}

/******** Flow control ********/

bool Mesodyn::mesodyn() {
  pepareForCalculations();
  if (success) {
    cout << "Mesodyn is all set, starting calculations.." << endl;
    for (int t = 0; t < timesteps; t++) {      // get segment potentials by iteration so that they match the given rho.
      New[0]->Solve(&rho[0], &dummyVector[0]); // get segment potentials by iteration so that they match the given rho.
      onsagerCoefficient();
      langevinFlux();
      updateDensity();
    }
  }
  return success;
}

void Mesodyn::pepareForCalculations() {
  componentNo = findComponentNo(); //find how many compontents there are (e.g. head, tail, solvent)
  size = findDensityVectorSize();  //find out how large the density vector is (needed for sizing the flux vector)
                                   //which will be 1 flux per lattice site per component per dimension
  buildRho();

  //find the index ranges where each component / dimension is located in the vector that
  //contains phi's and u's (which is a 3D lattice mapped onto 1D vector (Row-major))
  if (componentNo <= 1) {
    cout << "WARNING: Not enough components found for Mesodyn, aborting!" << endl;
    ;
    abort();
  } else if (componentNo > 1 && componentNo < 4) {
    setNeighborIndices(xNeighbors, yNeighbors, zNeighbors);
    setComponentStartIndices(component);
  } else {
    cout << "Unable to do Mesodyn for " << componentNo << " components, aborting!" << endl;
    ;
    abort();
  }
  //alocate memory for the fluxes of all components in all dimensions
  J.resize(size);
}

//defaults to homogeneous system for now
void Mesodyn::buildRho() {
  rho.reserve(size);
  for (int i = 0; i < size; ++i) {
    rho[i] = 0.5;
  }

  int steps = size / Lat[0]->MX;

  for (int i = 1; i < componentNo; ++i) {
    phi.push_back(&rho[steps * i]);
  }
}
void Mesodyn::abort() {
  //Once false is returned, Mesodyn automatically quits to main.
  success = false;
}

/****** Functions for handling indices in a 1D vector that contains a 3D lattice with multiple components *******/

int Mesodyn::findComponentNo() {
  component.reserve(In[0]->MolList.size());
  return In[0]->MolList.size();
}

int Mesodyn::findDensityVectorSize() {
  //  TODO: boundary conditions
  return componentNo * Lat[0]->MX * Lat[0]->MY * Lat[0]->MZ;
}

void Mesodyn::setNeighborIndices(vector<int>& xNeighbors, vector<int>& yNeighbors, vector<int>& zNeighbors) {
  //  TODO: boundary conditions
  zNeighbors[0] = -1;
  zNeighbors[1] = 1;

  yNeighbors[0] = -Lat[0]->MX;
  yNeighbors[1] = Lat[0]->MX;

  xNeighbors[0] = -Lat[0]->MX * Lat[0]->MY;
  xNeighbors[1] = Lat[0]->MX * Lat[0]->MY;
}

void Mesodyn::setComponentStartIndices(vector<int>& component) {
  component[0] = 0;
  component[1] = Lat[0]->MX * Lat[0]->MY * Lat[0]->MZ;
}

/******** Calculations ********/
void Mesodyn::onsagerCoefficient() {

  //TODO: maybe do this inline to preserve memory
  //TODO: optimize by doing all reserves once in Malloc.
  L.reserve(combinations(componentNo, 2) * Lat[0]->M);

  vector<Real>::iterator lIterator;
  lIterator = L.begin();

  //TODO: Untested code
  //all combinations of phi multiplications
  for (unsigned int i = 0; i < phi.size() - 1; ++i) {
    for (unsigned int j = i + 1; j < phi.size(); ++j) {
      //at every coordinate (pointer arithmatic)
      for (int xyz = 0; xyz < (Lat[0]->M); ++xyz) {
        *lIterator = (*(phi[i] + xyz) * *(phi[j] + xyz));
        lIterator++;
      }
    }
  }
}

void Mesodyn::langevinFlux() {
  vector<Real> u(size); //segment potential A

  //TODO: boundary condition in lattice?
  int z = 1;
  u[z] = New[0]->xx[z]; //which alphas? component a & b or one component at two sites?

  for (int i; i < componentNo; ++i) {
    J[z] = -D * L[z] + L[z + xNeighbors[0]];
  }

  /*for (int z = 1; z < size; z++) {
    gaussianNoise(dummyMean, dummyStdev, 2);
    //segment "chemical potential" gradient
    u[z] = alphaA[z] - alphaB[z];
    J[z] = (-D * (((L[z] + L[z + 1]) * (u[z + 1] - u[z])) - ((L[z - 1] + L[z]) * (u[z] - u[z - 1]))) + noise[0]);
    J[size + z] = (-D * (((L[z] + L[z + 1]) * (-u[z + 1] - -u[z])) - ((L[z - 1] + L[z]) * (-u[z] - -u[z - 1]))) + noise[1]);
  }*/
}

void Mesodyn::updateDensity() {
  //old density + langevinFluxTwo
}

/* Generates a vector of length count, contianing gaussian noise of given mean, standard deviation.
	 Noise is stored in vector<Real> Mesodyn::noise
	 Possible errors: What if count > sizeof(unsinged long)?
	 Called by langevinFlux()
*/
void Mesodyn::gaussianNoise(Real mean, Real stdev, unsigned long count) {

  random_device generator;

  seed_seq seed("something", "something else");

  //Mersenne Twister 19937 bit state size PRNG
  mt19937 prng(seed);

  normal_distribution<> dist(mean, stdev);

  this->noise.resize(count);

  for (unsigned int i = 0; i < count; ++i) {
    this->noise[i] = dist(prng);
  }

  /* Debugging code (output value in all elements):
for (auto const &element: mesodyn.thisNoise)
				std::cout << element << ' ';
*/
}

/******* Tools ********/
int Mesodyn::factorial(int n) {
  if (n > 1) {
    return n * factorial(n - 1);
  } else
    return 1;
}

int Mesodyn::combinations(int n, int k) {
  return factorial(n) / (factorial(n - k) * factorial(k));
}

void Mesodyn::PutParameter(string new_param) {
  KEYS.push_back(new_param);
}
string Mesodyn::GetValue(string parameter) {
  int i = 0;
  int length = PARAMETERS.size();
  while (i < length) {
    if (parameter == PARAMETERS[i]) {
      return VALUES[i];
    }
    i++;
  }
  return "";
}
void Mesodyn::push(string s, Real X) {
  Reals.push_back(s);
  Reals_value.push_back(X);
}
void Mesodyn::push(string s, int X) {
  ints.push_back(s);
  ints_value.push_back(X);
}
void Mesodyn::push(string s, bool X) {
  bools.push_back(s);
  bools_value.push_back(X);
}
void Mesodyn::push(string s, string X) {
  strings.push_back(s);
  strings_value.push_back(X);
}
void Mesodyn::PushOutput() {
  strings.clear();
  strings_value.clear();
  bools.clear();
  bools_value.clear();
  Reals.clear();
  Reals_value.clear();
  ints.clear();
  ints_value.clear();
}
Real* Mesodyn::GetPointer(string s) {
  //vector<string> sub;
  //nothing yet
  return NULL;
}
int Mesodyn::GetValue(string prop, int& int_result, Real& Real_result, string& string_result) {
  int i = 0;
  int length = ints.size();
  while (i < length) {
    if (prop == ints[i]) {
      int_result = ints_value[i];
      return 1;
    }
    i++;
  }
  i = 0;
  length = Reals.size();
  while (i < length) {
    if (prop == Reals[i]) {
      Real_result = Reals_value[i];
      return 2;
    }
    i++;
  }
  i = 0;
  length = bools.size();
  while (i < length) {
    if (prop == bools[i]) {
      if (bools_value[i])
        string_result = "true";
      else
        string_result = "false";
      return 3;
    }
    i++;
  }
  i = 0;
  length = strings.size();
  while (i < length) {
    if (prop == strings[i]) {
      string_result = strings_value[i];
      return 3;
    }
    i++;
  }
  return 0;
}
