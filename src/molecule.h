#ifndef MOLECULExH
#define MOLECULExH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "lattice.h"
#include "tools_host.h"
class Molecule {
public:
	Molecule(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::string);
virtual ~Molecule();

	std::string name;
	bool all_molecule;
	const Input* In;
	std::span<const std::unique_ptr<Segment>> Seg;
	Lattice* lat;
	std::vector<int> MolMonList;
	int start;
	Real KStiff;
	Real Mu;
	Real theta;
	Real phibulk;
	std::string freedom;
	MoleculeType MolType;
	Real n;
	Real GN;
	Real norm;
	int chainlength,N;
	std::vector<int> Gnr; //generation-number
	std::vector<int> first_s;
	std::vector<int> last_s;
	std::vector<int> first_b;
	std::vector<int> last_b;
	std::vector<int> mon_nr;
	std::vector<int> n_mon;
	std::vector<int> molmon_nr;
	std::vector<Real> mu_state;
	Real k_stiff;
	std::vector<Real> phi;
	std::vector<Real> phi_ranked;
	std::vector<Real> phitot;
	std::vector<Real> Gg_f;
	std::vector<Real> Gg_b;
	std::vector<Real> UNITY;
	std::vector<Real> P;
	int size;
    int Markov;
	Real B;
	ParameterStore OUTPUT;
	void PushOutput();
	std::span<Real> GetPointer(int);

	bool CheckInput(int,bool);
	bool ExpandBrackets(std::string&);
	bool Interpret(std::string,int);
	bool GenerateTree(std::string,int,int&,std::vector<int>,std::vector<int>);
	bool Decomposition(std::string);
	bool MakeMonList(void);
	bool IsPinned(void);
	bool IsCharged(void);
	Real Charge(void);
	bool HasOutputProperty(const std::string&) const;
	void DeAllocateMemory(void);
	void AllocateMemory(void);
	bool PrepareForCalculations(std::span<const Real>);
	virtual bool ComputePhi();
	virtual Real fraction(int);

	Real* propagate_forward(Real*,int&,int,int,int);
	void propagate_backward(Real*,int&,int,int,int);
	Real* propagate_forward(Real*,int&,int,Real*,int,int);
	void propagate_backward(Real*,int&,int,Real*,int&,int);

};

#endif
