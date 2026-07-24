#ifndef SAMPLE_CONFIGURATION_H
#define SAMPLE_CONFIGURATION_H

#include "input.h"
#include "molecule.h"

#include <cstdint>

class SampleConfiguration {
public:
	SampleConfiguration(std::string, const Molecule*, uint64_t, int);

	std::string name;
	std::string molecule_name;
	uint64_t seed;
	int size;
	ParameterStore xyz;

	// Generate an equilibrium layer path from the converged propagators, then lift it to representative XYZ coordinates.
	bool Generate();

private:
	const Molecule* molecule;
};

namespace sample_configuration {
bool Load(const Input&, Lattice*, std::span<const std::unique_ptr<Molecule>>, int, std::vector<SampleConfiguration>&);
}

#endif
