#ifndef MOL_CLAMPxH
#define MOL_CLAMPxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools_host.h"
class mol_clamp : public Molecule
{
	public: mol_clamp(const Input*,Lattice*,vector<Segment*>,string);
	~mol_clamp();

	bool ComputePhi();
};

#endif
