#ifndef OUTPUTxH
#define OUTPUTxH
#include "namics.h"
#include "input.h"
#include "lattice.h"
#include "segment.h"
#include "state.h"
#include "reaction.h"
#include "molecule.h"
#include "system.h"
#include "solve_scf.h"
#include "json_writer.h"

class Output {
public:
	Output(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::span<const std::unique_ptr<State>>,std::span<const std::unique_ptr<Reaction>>,std::span<const std::unique_ptr<Molecule>>,System*,Solve_scf*,std::string);

~Output();

	std::string name;
	const Input* In;
	Lattice* lat;
	std::span<const std::unique_ptr<Segment>> Seg;
	std::span<const std::unique_ptr<State>> Sta;
	std::span<const std::unique_ptr<Reaction>> Rea;
	std::span<const std::unique_ptr<Molecule>> Mol;
	System* Sys;
	Solve_scf* New;
	std::shared_ptr<io::json::JsonWriter> json_writer;
	int start;
	bool write_bounds;
	bool append;
	bool write;
	std::string sep;

  	std::vector<std::string> ints;
 	std::vector<std::string> Reals;
 	std::vector<std::string> bools;
	std::vector<std::string> strings;
  	std::vector<Real> Reals_value;
  	std::vector<int> ints_value;
  	std::vector<bool> bools_value;
  	std::vector<std::string> strings_value;

	std::vector<std::string> OUT_key;
	std::vector<std::string> OUT_name;
	std::vector<std::string> OUT_prop;


	std::vector<std::string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutParameter(std::string);
	std::string GetValue(std::string);
	bool Load();
	void WriteOutput(int);
	int GetValue(std::string, std::string, std::string, int&, Real&, std::string&);
	std::span<Real> GetPointer(std::string, std::string, std::string);
	std::span<int> GetPointerInt(std::string, std::string, std::string);
	int GetValue(std::string, std::string, int& , Real& , std::string&);
	void push(std::string, Real);
	void push(std::string, int);
	void push(std::string, bool);
	void push(std::string, std::string);

};
#endif
