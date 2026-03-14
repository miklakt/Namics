#ifndef REACTIONxH
#define REACTIONxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "state.h"

#include <cmath>

class Reaction {
public:
	Reaction(const Input*,std::span<const std::unique_ptr<Segment>>,std::span<const std::unique_ptr<State>>,std::string);

	~Reaction();
	std::span<const std::unique_ptr<State>> Sta;
	std::span<const std::unique_ptr<Segment>> Seg;
	std::vector<int> Sto;
	std::vector<int> State_nr;
	std::vector<int> State_in_seg_nr; 
	std::vector<int> Seg_nr; 
	Real K;
	Real pK;
	std::string equation;
	
	std::string name; 
	const Input* In;
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
	bool CheckInput(int);
	void PutParameter(std::string); 
	std::string GetValue(std::string); 
	Real ChemIntBulk(const State&);
	Real pKeff();
	Real Residual_value();
	bool GuessAlpha();
	bool PutAlpha(Real);
};
#endif
