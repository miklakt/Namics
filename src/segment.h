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
	Lattice* Lat;
	Lattice* lat;
	vector<int> constraint_z;
	vector<Real> constraint_phi;
	vector<Real> constraint_beta;
	bool constraints;

	vector<string> chi_name;
	vector<Real> chi;
	int n_seg;
	int seg_nr;
	bool unique;
	int seg_nr_of_copy;
	int state_nr_of_copy;
	bool prepared;

	Real theta_exc;
	Real M1,M2,Fl;
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

	string filename;
	string copy_of;
	string s_freedom;
	bool block;
	bool all_segment;
	int ns;

	int n_pos;
	int* r;
	int start;
	int var_pos;
	int frozen_at_bound;
	int used_in_mol_nr;
	Real phi_LB_X;
	Real phi_UB_X;
	Real phi_LB_Y;
	Real phi_UB_Y;

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
	Real Get_g(int ) ;
	void Put_beta(int, Real );
	void PutContraintBC ();


	string GetOriginal();

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
	Real B;
	Real J;

	std::vector<string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutChiKEY(string);
	string GetValue(string);
	string GetFreedom();
	Real PinnedVolume();
	bool IsFree();
	bool IsPinned();
	bool IsFrozen();
	Real* GetMASK();
	Real* GetPhi();
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
