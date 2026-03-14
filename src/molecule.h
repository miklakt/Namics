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
	Lattice* lat;
	vector<int> MolMonList;
	int start;
	Real KStiff;
	Real Mu;
	Real theta;
	Real phibulk;
	string freedom;
	MoleculeType MolType;
	Real n;
	Real GN;
	Real norm;
	int chainlength,N;
	vector<int> Gnr; //generation-number
	vector<int> first_s;
	vector<int> last_s;
	vector<int> first_b;
	vector<int> last_b;
	vector<int> mon_nr;
	vector<int> n_mon;
	vector<int> molmon_nr;
	vector<Real> mu_state;
	Real k_stiff;
	Real *phi;
	Real *rho;
	Real *H_phi;
	Real *phitot;
	Real *H_phitot;
	Real *Gg_f;
	Real *Gg_b;
	Real *UNITY;
	Real *P;
	int size;
    int Markov;
	Real B;

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

	std::vector<string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int,bool);
	void PutParameter(string);
	bool ExpandBrackets(string&);
	bool Interpret(string,int);
	bool GenerateTree(string,int,int&,vector<int>,vector<int>);
	bool Decomposition(string);
	string GetValue(string);
	bool MakeMonList(void);
	bool IsPinned(void);
	bool IsCharged(void);
	Real Charge(void);
	void DeAllocateMemory(void);
	void AllocateMemory(void);
	bool PrepareForCalculations(Real*);
	virtual bool ComputePhi();
	virtual Real fraction(int);

	Real* propagate_forward(Real*,int&,int,int,int);
	void propagate_backward(Real*,int&,int,int,int);
	Real* propagate_forward(Real*,int&,int,Real*,int,int);
	void propagate_backward(Real*,int&,int,Real*,int&,int);

};

#endif
