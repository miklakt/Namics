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
#include "sample_configuration.h"
#include "json_writer.h"

class Output {
public:
	Output(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::span<const std::unique_ptr<State>>,std::span<const std::unique_ptr<Reaction>>,std::span<const std::unique_ptr<Molecule>>,std::span<const SampleConfiguration>,System*,Solve_scf*,std::string);

	std::string name;
	const Input* In;
	Lattice* lat;
	std::span<const std::unique_ptr<Segment>> Seg;
	std::span<const std::unique_ptr<State>> Sta;
	std::span<const std::unique_ptr<Reaction>> Rea;
	std::span<const std::unique_ptr<Molecule>> Mol;
	std::span<const SampleConfiguration> Samples;
	System* Sys;
	Solve_scf* New;
	io::json::JsonWriter json_writer;
	int start;
	bool write_bounds;
	bool append;
	bool write;
	std::string sep;
	ParameterStore items;

	bool CheckInput(int);
	bool Load();
	void WriteOutput(int);
};
#endif
