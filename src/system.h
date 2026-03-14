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
	System(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::span<const std::unique_ptr<State>>,std::span<const std::unique_ptr<Reaction>>,std::span<const std::unique_ptr<Molecule>>,std::string);

~System();

	std::string name;
	bool all_system;
	const Input* In;
	std::vector<Real> CHI;
	std::span<const std::unique_ptr<Segment>> Seg;
	std::span<const std::unique_ptr<State>> Sta;
	std::span<const std::unique_ptr<Reaction>> Rea;
	std::span<const std::unique_ptr<Molecule>> Mol;
	Lattice* lat;
	std::vector<int> SysMonList;
	std::vector<int> ItMonList;
	std::vector<int> ItStateList;
	std::vector<int> FrozenList;
	std::vector<int> XmolList;
	std::vector<int> XstateList_1;
	std::vector<int> XstateList_2;
	std::vector<int> Xn_1;
	Real FreeEnergy;
	Real GrandPotential;

	std::vector<Real> phitot;
	std::vector<Real> KSAM;
	std::vector<Real> eps;
	std::vector<Real> psi;
	std::vector<Real> EE;
	std::vector<Real> E;
	std::vector<Real> psiMask;
	bool fixedPsi0;
	bool grad_epsilon;
	std::vector<Real> q;
	std::vector<Real> GrandPotentialDensity;
	std::vector<Real> FreeEnergyDensity;
	std::vector<Real> alpha;
	std::vector<Real> TEMP;
	bool first_pass;
	int n_mol;
	int solvent;
	int neutralizer;
	bool charged;
	std::string initial_guess;
	std::string guess_inputfile;
	bool write_initial_guess;
	int start;

	std::vector<std::string> ints;
	std::vector<std::string> Reals;
	std::vector<std::string> bools;
	std::vector<std::string> strings;
	std::vector<Real> Reals_value;
	std::vector<int> ints_value;
	std::vector<bool> bools_value;
	std::vector<std::string> strings_value;
	void push(std::string,Real);
	void push(std::string,int);
	void push(std::string,bool);
	void push(std::string,std::string);
	void PushOutput();
	std::span<Real> GetPointer(std::string);
	std::span<int> GetPointerInt(std::string);
	int GetValue(std::string,int&,Real&,std::string&);

	std::vector<std::string> KEYS;
	ParameterStore PARAMETERS;
	void DeAllocateMemory();
	bool CheckInput(int);
	void PutParameter(std::string);
	std::string GetValue(std::string);
	bool CheckChi_values(int);
	bool MakeItsLists();
	bool IsUnique(int,int);
	void AllocateMemory();
	bool PrepareForCalculations(bool);
	bool generate_mask();
	bool ComputePhis(Real);
	void ComputePhis(std::span<const Real>,bool,Real);
	bool PutU(std::span<const Real>);
	void Classical_residual(std::span<const Real>,std::span<Real>,Real,int);

	void DoElectrostatics(std::span<Real>,std::span<const Real>);
	bool CheckResults(bool);
	Real GetE(int,int);
	Real GetFreeEnergy();
	Real GetGrandPotential();
	bool CreateMu(int);
};
#endif
