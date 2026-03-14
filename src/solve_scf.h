#ifndef SOLVE_SCFxH
#define SOLVE_SCFxH
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
#include <Eigen/Core>
#include <LBFGS.h>
typedef Eigen::Matrix<Real,Eigen::Dynamic,1> Vector;
typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> Matrix;
using namespace LBFGSpp;


class Solve_scf : public SFNewton {
public:
	Solve_scf() {};

	Solve_scf(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::span<const std::unique_ptr<State>>,std::span<const std::unique_ptr<Reaction>>,std::span<const std::unique_ptr<Molecule>>,System*,std::string);

	~Solve_scf();

	std::string name;
	const Input* In;
	System* Sys;
	std::span<const std::unique_ptr<Segment>> Seg;
	Lattice* lat;
	std::span<const std::unique_ptr<Molecule>> Mol;
	std::span<const std::unique_ptr<State>> Sta;
	std::span<const std::unique_ptr<Reaction>> Rea;

	int start;
	std::string SCF_method;
	std::string stop_criterion;
	int iv;
	int m, restart_DIIS;
	bool all;

	Real tolerance;
	Real deltamax,deltamin;

	int iterationlimit;

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
	int GetValue(std::string,int&,Real&,std::string&);
	enum iteration_method {HESSIAN,PSEUDOHESSIAN,diis,LBFGS};
	enum inner_iteration_method {super,proceed};
	enum gradient_method {classical, WEAK};
	iteration_method solver;
	gradient_method gradient;
	inner_iteration_method control;


	std::vector<Real> xx;
	std::vector<Real> yy;
	std::vector<int> SIGN;

	std::vector<std::string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutParameter(std::string);
	std::string GetValue(std::string);
	void Copy(std::span<Real>,std::span<const Real>,int,int,int,int);
	bool Guess(std::span<const Real>,std::string,std::vector<std::string>,std::vector<std::string>,bool,int,int,int,int);

	bool Solve(bool);

	void DeAllocateMemory();
	void AllocateMemory();
	void residuals(Real*,Real*);

	void inneriteration(Real*,Real*,Real*,Real,Real&,Real,int);

};
#endif
