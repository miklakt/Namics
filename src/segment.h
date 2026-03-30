#ifndef SEGMENTxH
#define SEGMENTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "tools.h"
class Segment {
public:
	Segment(const Input*,Lattice*,std::string,int,int);

~Segment();

	std::string name;
	const Input* In;
	Lattice* lat;

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

	bool all_segment;
	int ns;

	int start;
	int var_pos;
	int frozen_at_bound;
	ParameterStore OUTPUT;
	void PushOutput();
	std::span<Real> GetPointer(int);
	bool LoadExternalPotential();

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

	bool CheckInput(int);
	Real PinnedVolume();
	void DeAllocateMemory();
	void AllocateMemory();
	bool PrepareForCalculations(std::span<const Real>,bool);
	bool ParseFreedoms();
	void UpdateValence(Real*,std::span<Real>,std::span<Real>,std::span<Real>,bool);
	int AddState(int,Real,Real,bool);
	void SetPhiSide();
	bool PutAlpha(Real);
	bool CanBeReached(int, int, int, int);
};
#endif
