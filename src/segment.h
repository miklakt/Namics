#ifndef SEGMENTxH
#define SEGMENTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "tools_host.h"
class Segment {
public:
	Segment(const Input*,Lattice*,string,int,int);

~Segment();

	string name;
	const Input* In;
	Lattice* lat;

	vector<string> chi_name;
	vector<Real> chi;
	int n_seg;
	int seg_nr;
	bool unique;
	int seg_nr_of_copy;
	int state_nr_of_copy;

	Real theta_exc;
	Real epsilon;
	Real valence;
	Real PSI0;
	bool fixedPsi0;
	Real phibulk;
	string freedom;
	vector<int>state_change;
	vector<Real>state_valence;
	vector<int>state_id;
	vector<string>state_name;
	vector<int>state_nr;
	vector<Real>state_alphabulk;
	vector<Real>state_phibulk;
	vector<Real>state_theta;

	bool block;
	bool all_segment;
	int ns;

	int n_pos;
	std::array<int, 6> r;
	int start;
	int var_pos;
	int frozen_at_bound;

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
	bool LoadExternalPotential();

	int* H_P;
	Real* H_MASK;
	Real* H_u;
	Real* H_phi;
	Real* H_u_ext;

	Real* H_phi_state;
	Real* H_alpha;
	Real* H_ALPHA;

	int* P;
	Real* MASK;
	Real* G1;
	Real* phi;
	Real* phi_state;
	Real* phi_side;
	Real* u;
	Real* u_ext;

	Real* alpha;//fraction of segment in specfied state
	Real* ALPHA; //Lagrange parameter per segement for steady state
	int ItState;

	std::vector<string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutChiKEY(string);
	string GetValue(string);
	Real PinnedVolume();
	void DeAllocateMemory();
	void AllocateMemory();
	bool PrepareForCalculations(Real*,bool);
	bool ParseFreedoms(bool&);
	void UpdateValence(Real*,Real*,Real*,Real*,bool);
	int AddState(int,Real,Real,bool);
	void SetPhiSide();
	bool PutAlpha(Real);
	bool CanBeReached(int, int, int, int);
};
#endif
