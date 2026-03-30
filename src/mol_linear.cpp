#include "molecule.h"
#include "mol_linear.h"


mol_linear::mol_linear(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_, std::string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_linear::~mol_linear() {
}



bool mol_linear::ComputePhi() {
NAMICS_DBG("ComputePhi in mol_linear " << std::endl);
	int b0 = first_b[0];
	int bN = last_b[0];
	int M=lat->M;
	bool success=true;
	int s=0;
	Real* Glast=NULL;
	for (int b = b0; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1.data(),s,b,0,M);

	GN=lat->ComputeGN(Glast,M);

	s--;
	for (int b = bN ; b >= b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1.data(),s,b,M);
	return success;
}
