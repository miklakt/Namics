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
	ParameterStore OUTPUT;
	void PushOutput();

	bool CheckInput(int);
	Real ChemIntBulk(const State&);
	Real pKeff();
	Real Residual_value();
	bool GuessAlpha();
	bool PutAlpha(Real);
};
#endif
