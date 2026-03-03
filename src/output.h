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
#include "io_utils.h"
#include "json_writer.h"
#include <climits>
#include <unistd.h>

class Output {
public:
	Output(const Input*,Lattice*,vector<Segment*>,vector<State*>,vector<Reaction*>,vector<Molecule*>,System*,Solve_scf*,string,int,int);

~Output();

	string name;
	const Input* In;
	Lattice* Lat;
	Lattice* lat;
	vector<Segment*> Seg;
	vector<State*> Sta;
	vector<Reaction*> Rea;
	vector<Molecule*> Mol;
	System* Sys;
	Solve_scf* New;
	std::shared_ptr<io::legacy::Writer> writer;
	std::shared_ptr<io::json::JsonWriter> json_writer;
	int n_output;
	int subl;
	int n_starts;
	int start;
	int output_nr;
	bool write_bounds;
	bool append;
	bool write;
	bool input_error;
	bool DOS;
	int first;
	string output_folder;
	string bin_folder;
	bool use_output_folder;
	string write_option;
	string sep;

  	vector<string> ints;
 	vector<string> Reals;
 	vector<string> bools;
 	vector<string> strings;
  	vector<Real> Reals_value;
	vector<int > SizeVectorReal;
	vector<int > SizeVectorInt;
	vector<Real*> PointerVectorReal;
	vector<int*> PointerVectorInt;
	vector<int> pointer_size;
  	vector<int> ints_value;
  	vector<bool> bools_value;
  	vector<string> strings_value;

	std::vector<string> OUT_key;
	std::vector<string> OUT_name;
	std::vector<string> OUT_prop;


	std::vector<string> KEYS;
	ParameterStore PARAMETERS;
	bool CheckInput(int);
	void PutParameter(string);
	string GetValue(string);
	bool Load();
	void WriteOutput(int);
	int GetValue(string, string, string, int&, Real&, string&);
	Real* GetPointer(string, string, string, int&);
	int* GetPointerInt(string, string, string, int&);
	int GetValue(string, string, int& , Real& , string&);
	void push(string, Real);
	void push(string, int);
	void push(string, bool);
	void push(string, string);

};
#endif
