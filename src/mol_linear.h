#ifndef MOL_LINEARxH
#define MOL_LINEARxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "lattice.h"
#include "tools_host.h"
class mol_linear : public Molecule
{
	public: mol_linear(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::string);
	~mol_linear();

	bool ComputePhi();

};

#endif
