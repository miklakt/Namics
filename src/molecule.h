#ifndef MOLECULExH
#define MOLECULExH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "lattice.h"
#include "tools_host.h"
class Molecule {
public:
	Molecule(const Input*,Lattice*,vector<Segment*>,string);
virtual ~Molecule();

	string name;
	bool all_molecule;
	const Input* In;
	vector<Segment*> Seg;
	Lattice* Lat;
	Lattice* lat;
	vector<int> MolMonList;
	vector<int> FillRangesList;
	int start;
	int n_mol;
	int mol_nr;
	Real KStiff;
	Real Mu;
	Real theta;
	Real R_Gibbs;
	Real theta_Gibbs;
	Real theta_range,n_range;
	Real phibulk;
	string freedom;
	MoleculeType MolType;
	Real n;
	Real GN,GN1,GN2;
	Real norm;
	Real phi1,phiM,width,Dphi,pos_interface,phi_av;
	int chainlength,N;
	bool save_memory;
	string composition;
	vector<int> Gnr; //generation-number
	vector<int> first_s;
	vector<int> last_s;
	vector<int> first_b;
	vector<int> last_b;
	vector<int> mon_nr;
	vector<int> n_mon;
	vector<int> molmon_nr;
	vector<int> memory;
	vector<int> last_stored;
	vector<Real> mu_state;
	vector<Real> block;
	Real k_stiff;
	Real *phi;
	Real *rho;
	Real *H_phi;
	Real *R_mask;
	Real *phitot;
	Real *H_phitot;
	Real *Gg_f;
	Real *Gg_b;
	Real *Gs;
	Real *UNITY;
	Real *P;
	int size;
    int Markov;
	bool ring;
	Real J;
	Real B;
	Real Delta_MU;

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
	void PutTheta(Real);
	bool ComputeWidth();
	Real ComputeGibbs(Real);
	void SetThetaBlocks(int);
	void NormPerBlock(int);
	bool Filling;

	std::vector<string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int,bool);
	void PutParameter(string);
	bool ExpandBrackets(string&);
	bool Interpret(string,int);
	bool GenerateTree(string,int,int&,vector<int>,vector<int>);
	bool Decomposition(string);
	int GetChainlength(void);
	Real Theta(void);
	string GetValue(string);
	int GetMonNr(string);
	bool MakeMonList(void);
	bool IsPinned(void);
	int GetPinnedSeg(void);
	bool IsCharged(void);
	Real Charge(void);
	void DeAllocateMemory(void);
	void AllocateMemory(void);
	bool PrepareForCalculations(Real*);
	bool ComputePhi(Real*,int);
	virtual bool ComputePhi();
	virtual Real fraction(int);

	Real* propagate_forward(Real*,int&,int,int,int);
	void propagate_backward(Real*,int&,int,int,int);
	Real* propagate_forward(Real*,int&,int,Real*,int,int);
	void propagate_backward(Real*,int&,int,Real*,int&,int);

};

#endif
