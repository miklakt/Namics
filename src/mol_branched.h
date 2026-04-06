#ifndef MOL_BRANCHEDxH
#define MOL_BRANCHEDxH
#include "molecule.h"

class mol_branched : public Molecule {
public:
	mol_branched(const Input*,Lattice*,std::span<const std::unique_ptr<Segment>>,std::string);
	~mol_branched() override = default;
	bool ComputePhi() override;

private:
	Real* ForwardBranch(int generation, int &s);
	void BackwardBranch(int generation, int &s);
};

#endif
