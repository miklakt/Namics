#ifndef MOL_BRANCHEDxH
#define MOL_BRANCHEDxH
#include "namics.h"
#include "input.h"
#include "segment.h"
#include "lattice.h"
#include "tools.h"
class mol_branched : public Molecule
{
	public: mol_branched(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::string);
	~mol_branched();

	Real* ForwardBra(int generation, int &s);
	void BackwardBra(int generation, int &s);




	bool ComputePhi();

};

#endif
