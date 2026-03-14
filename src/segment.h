#ifndef SEGMENTxH
#define SEGMENTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "tools_host.h"
class Segment {
public:
	Segment(const Input*,Lattice*,std::string,int,int);

~Segment();

	std::string name;
	const Input* In;
	Lattice* lat;

	std::vector<std::string> chi_name;
	std::vector<Real> chi;
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
	std::string freedom;
	std::vector<int>state_change;
	std::vector<Real>state_valence;
	std::vector<int>state_id;
	std::vector<std::string>state_name;
	std::vector<int>state_nr;
	std::vector<Real>state_alphabulk;
	std::vector<Real>state_phibulk;
	std::vector<Real>state_theta;

	bool block;
	bool all_segment;
	int ns;

	int n_pos;
	std::array<int, 6> r;
	int start;
	int var_pos;
	int frozen_at_bound;

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
	bool LoadExternalPotential();

	std::vector<int> P;
	std::vector<Real> MASK;
	std::vector<Real> G1;
	std::vector<Real> phi;
	std::vector<Real> phi_state;
	std::vector<Real> phi_side;
	std::vector<Real> u;
	std::vector<Real> u_ext;

	std::vector<Real> alpha;//fraction of segment in specfied state
	std::vector<Real> ALPHA; //Lagrange parameter per segement for steady state
	int ItState;

	std::vector<std::string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutChiKEY(std::string);
	std::string GetValue(std::string);
	Real PinnedVolume();
	void DeAllocateMemory();
	void AllocateMemory();
	bool PrepareForCalculations(std::span<const Real>,bool);
	bool ParseFreedoms(bool&);
	void UpdateValence(Real*,std::span<Real>,std::span<Real>,std::span<Real>,bool);
	int AddState(int,Real,Real,bool);
	void SetPhiSide();
	bool PutAlpha(Real);
	bool CanBeReached(int, int, int, int);
};
#endif
