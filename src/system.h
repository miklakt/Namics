#ifndef SYSTEMxH
#define SYSTEMxH
#include "namics.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "molecule.h"
#include "lattice.h"

class System {
public:
	System(const Input*,Lattice*,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,string);

~System();

	string name;
	bool all_system;
	const Input* In;
	Real* CHI;
	vector<Segment*> Seg;
	vector<State*> Sta;
	vector<Reaction*> Rea;
	vector<Molecule*> Mol;
	Lattice* Lat;
	Lattice* lat;
	vector<int> SysMonList;
	vector<int> ItMonList;
	vector<int> ItStateList;
	vector<int> FrozenList;
	vector<int> SysTagList;
	vector<int> FillList;
	vector<int> XmolList;
	vector<int> XstateList_1;
	vector<int> XstateList_2;
	vector<int> DeltaMolList;
	vector<int> px;
	vector<int> py;
	vector<int> pz;
	vector<int> Xn_1;
	Real FreeEnergy;
	Real GrandPotential;
	Real pos_interface;

	Real* phitot;
	Real* KSAM;
	Real* FILL;
	Real* eps;
	Real* H_psi;
	Real* H_q;
	Real* psi;
	Real* EE;
	Real* E;
	Real* psiMask;
	bool fixedPsi0;
	bool grad_epsilon;
	bool constraintfields;
	Real* q;
	Real* H_GrandPotentialDensity;
	Real* H_FreeEnergyDensity;
	Real* H_alpha;
	Real* GrandPotentialDensity;
	Real* FreeEnergyDensity;
	Real* alpha;
	Real* TEMP;
	Real* H_BETA;
	Real* BETA;
	Real phi_ratio;
	Real* H_beta;
	Real* beta;
	bool first_pass;
	int n_mol;
	int solvent;
	int neutralizer;
	int tag_segment;
	Real volume;
	bool charged;
	bool local_solution;
	bool do_blocks;
	int split;
	string initial_guess;
	string guess_inputfile;
	bool write_initial_guess;
	string ConstraintType;
	string delta_inputfile;
	int start;
	int extra_constraints;
	Real old_residual;
	int progress;


	bool prepared;
	bool Filling;

	vector<string> ints;
	vector<string> Reals;
	vector<string> bools;
	vector<string> strings;
	vector<Real> Reals_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;
	void push(string,Real);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	Real* GetPointer(string,int&);
	int* GetPointerInt(string,int&);
	int GetValue(string,int&,Real&,string&);
	string CalculationType; // {micro_emulsion,micro_phasesegregation};

	std::vector<string> KEYS;
	ParameterStore PARAMETERS;
	void DeAllocateMemory();
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);
	string GetMonName(int );
	int GetMonNr(string);
	bool CheckChi_values(int);
	bool MakeItsLists();
	bool IsCharged();
	bool IsUnique(int,int);
	void AllocateMemory();
	bool PrepareForCalculations(bool);
	bool generate_mask();
	bool ComputePhis(Real);
	void ComputePhis(Real*,bool,Real);
	bool Put_U(Real*);
	bool PutU(Real*);
	void Classical_residual(Real* ,Real*,Real,int,int);

	void DoElectrostatics(Real*,Real*);
	bool CheckResults(bool);
	Real GetE(int,int);
	Real GetFreeEnergy();
	Real GetGrandPotential();
	Real GetSpontaneousCurvature();
	Real GetKBar();
	bool CreateMu(int);
};
#endif
