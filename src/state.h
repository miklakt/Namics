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
	std::vector<std::string> chi_name;
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
	void PutChiKEY(std::string);
	std::vector<std::string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutParameter(std::string); 
	std::string GetValue(std::string); 
};
#endif
