#ifndef SOLVE_SCFxH
#define SOLVE_SCFxH
#include <limits>
#include "namics.h"
#include "input.h"
#include "system.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "lattice.h"
#include "molecule.h"
#include "tools_host.h"
#include "sfnewton.h"
#include <functional>
#include <Eigen/Core>
#include "LBFGS.h"
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> Vector;
typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> Matrix;
using namespace LBFGSpp;


class Solve_scf : public SFNewton {
public:
	Solve_scf() {};

	Solve_scf(const Input*,Lattice*,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,System*,string);

	~Solve_scf();

	enum rescue {
		NONE,
		ZERO,
		M,
		DELTA_MAX,
	};

	bool attempt_DIIS_rescue();
	rescue rescue_status;

	string name;
	const Input* In;
	System* Sys;
	vector<Segment*> Seg;
	Lattice* Lat;
	Lattice* lat;
	vector<Molecule*> Mol;
	vector<State*> Sta;
	vector<Reaction*> Rea;

	int start;
	string SCF_method;
	string gradients;
	string StoreFileGuess;
	string ReadFileGuess;
	string stop_criterion;
	int iv;
	int m, restart_DIIS;
	SolverType SolType;
	bool all;

	Real tolerance;
	Real deltamax,deltamin;

	int iterationlimit;
	bool value_e_info;
	bool value_s_info;
	int value_i_info;
	Real* temp_alpha;

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
	int GetValue(string,int&,Real&,string&);
	enum iteration_method {HESSIAN,PSEUDOHESSIAN,PICARD,diis,conjugate_gradient,LBFGS,BRR};
	enum inner_iteration_method {super,proceed};
	enum gradient_method {classical, Picard, WEAK};
	iteration_method solver;
	gradient_method gradient;
	inner_iteration_method control;


	Real *xx;
	Real *yy;
	vector<Real> x_storage;
	int *SIGN;
	Real* alpha;

	std::vector<string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);
	void Copy(Real*,Real*,int,int,int,int);
	bool Guess(Real*,string,vector<string>,vector<string>,bool,int,int,int,int);

	bool Solve(bool);

	void DeAllocateMemory();
	void AllocateMemory();
	bool PrepareForCalculations(void);
	void ComputePhis(bool);
	bool PutU();
	bool Put_U();
	void residuals(Real*,Real*);
	void gradient_log(Real*, int, int, int, int);
	void gradient_quotient(Real*, int, int, int, int);
	void gradient_minus(Real*, int, int, int, int);
	function<void(Real*, int, int, int, int)> target_function;

	void inneriteration(Real*,Real*,Real*,Real,Real&,Real,int);

};
#endif
