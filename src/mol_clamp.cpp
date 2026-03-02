#include "molecule.h"
#include "mol_clamp.h"


mol_clamp::mol_clamp(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {}


mol_clamp::~mol_clamp() {
}


bool mol_clamp::ComputePhi(){
	NAMICS_DBG_THIS("ComputePhi for mol_clamp " + name << endl);	bool success=true;
	int M=lat->M;
	int m=0;
	if (freedom=="clamped") m=lat->m[Seg[mon_nr[0]]->clamp_nr];
	int blocks=mon_nr.size();
	std::fill_n(rho, m*n_box*MolMonList.size(), 0);
	int s=1;
	if (save_memory) {
		std::copy_n(mask1, m*n_box, Gs);
	} else {
		std::copy_n(mask1, m*n_box, Gg_f);
	}
	for (int i=1; i<blocks-1; i++) {
		lat->DistributeG1(std::span<const Real>(Seg[mon_nr[i]]->G1, static_cast<size_t>(M)),
		                 std::span<Real>(g1, static_cast<size_t>(m*n_box)),
		                 std::span<const int>(Bx, static_cast<size_t>(n_box)),
		                 std::span<const int>(By, static_cast<size_t>(n_box)),
		                 std::span<const int>(Bz, static_cast<size_t>(n_box)),
		                 n_box);

		propagate_forward(g1,s,i,0,m*n_box);
	}
	if (save_memory) {
		int k=last_stored[blocks-2];
		int N=memory[n_mon.size()-1];
		lat->propagate(Gg_f,mask2,k,N-1,m*n_box);
		lat->ComputeGN(std::span<Real>(gn, static_cast<size_t>(n_box)),
		               std::span<const Real>(Gg_f, static_cast<size_t>(n_box*m*N)),
		               std::span<const int>(H_Bx, static_cast<size_t>(n_box)),
		               std::span<const int>(H_By, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Bz, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Px2, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Py2, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Pz2, static_cast<size_t>(n_box)),
		               N-1,
		               n_box);
	} else {
		lat->propagate(Gg_f,mask2,s-1,s,m*n_box);
		lat->ComputeGN(std::span<Real>(gn, static_cast<size_t>(n_box)),
		               std::span<const Real>(Gg_f, static_cast<size_t>(n_box*m*chainlength)),
		               std::span<const int>(H_Bx, static_cast<size_t>(n_box)),
		               std::span<const int>(H_By, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Bz, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Px2, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Py2, static_cast<size_t>(n_box)),
		               std::span<const int>(H_Pz2, static_cast<size_t>(n_box)),
		               chainlength-1,
		               n_box);
	}
	s=chainlength-1;
	std::copy_n(mask2, m*n_box, Gg_b+(s%2)*m*n_box);
	if (save_memory) std::copy_n(Gg_b+(s%2)*m*n_box, m*n_box, Gg_b+((s-1)%2)*m*n_box);
	s--;
	for (int i=blocks-2; i>0; i--) {
		lat->DistributeG1(std::span<const Real>(Seg[mon_nr[i]]->G1, static_cast<size_t>(M)),
		                 std::span<Real>(g1, static_cast<size_t>(m*n_box)),
		                 std::span<const int>(Bx, static_cast<size_t>(n_box)),
		                 std::span<const int>(By, static_cast<size_t>(n_box)),
		                 std::span<const int>(Bz, static_cast<size_t>(n_box)),
		                 n_box);

		propagate_backward(g1,s,i,0,m*n_box);
	}
	for (size_t i = 1; i < MolMonList.size(); i++ )
	{
		lat->CollectPhi(std::span<Real>(phi+M*i, static_cast<size_t>(M)),
		               std::span<const Real>(gn, static_cast<size_t>(n_box)),
		               std::span<const Real>(rho+m*n_box*i, static_cast<size_t>(m*n_box)),
		               std::span<const int>(Bx, static_cast<size_t>(n_box)),
		               std::span<const int>(By, static_cast<size_t>(n_box)),
		               std::span<const int>(Bz, static_cast<size_t>(n_box)),
		               n_box);
	}

	return success;
}




