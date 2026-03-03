#ifndef MOL_LINEARxH
#define MOL_LINEARxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "lattice.h"
#include "tools_host.h"
class mol_linear : public Molecule
{
	public: mol_linear(const Input*,Lattice*,vector<Segment*>,string);
	~mol_linear();

	bool ComputePhi();

};

#endif
