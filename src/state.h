#ifndef STATExH
#define STATExH
#include "namics.h"
#include "input.h"
#include "segment.h"
class State {
public:
	State(const Input*,std::span<const std::unique_ptr<Segment>>,std::string);

	~State();
	const Input* In;
	std::span<const std::unique_ptr<Segment>> Seg;
	std::string name; 
	std::vector<Real> chi;
 
	bool unique;
	int seg_nr_of_copy;
	int state_nr_of_copy; 
	Real valence;
	std::string mon_name;
	int mon_nr; //number of the monomer to which the state belongs
	int state_nr; //number of the state of the monomer it belongs to
	int state_id; //state number in In->StateList
	std::string state_name;
	Real alphabulk;
	bool fixed;
	bool in_reaction;
	ParameterStore OUTPUT;
	void PushOutput();
	std::span<Real> GetPointer(int);
	bool CheckInput(int);
};
#endif
