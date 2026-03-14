#include "molecule.h"
#include "mol_linear.h"


mol_linear::mol_linear(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_linear::~mol_linear() {
}



bool mol_linear::ComputePhi() {
NAMICS_DBG("ComputePhi in mol_linear " << endl);
	int b0 = first_b[0];
	int bN = last_b[0];
	int M=lat->M;
	bool success=true;
	int unity=0;
	int s=0;
	Real* Glast=NULL;
	if (Markov ==2)
		for (int b = b0; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1,s,b,P,0,M);
	else
		for (int b = b0; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1,s,b,0,M);

	GN=lat->ComputeGN(Glast,Markov,M);

	s--;
	if (Markov==2)
		for (int b = bN ; b >= b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1,s,b,P,unity,M);
	else
		for (int b = bN ; b >= b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1,s,b,0,M);
	return success;
}


