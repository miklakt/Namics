#ifndef MOL_DENDxH
#define MOL_DENDxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "alias.h"
#include "lattice.h"
#include "tools_host.h"
class mol_dend : public Molecule
{
	public: mol_dend(const Input*,Lattice*,vector<Segment*>,string);
	~mol_dend();

	Real fraction(int);
	bool ComputePhi();
	bool BackAndForth();
	bool BackAndForth2ndO();
};

#endif
