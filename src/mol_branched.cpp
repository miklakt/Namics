#include "mol_branched.h"


mol_branched::mol_branched(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_, std::string name_) : Molecule(In_,Lat_,Seg_,name_) {

}

void mol_branched::BackwardBranch(int generation, int &s){
NAMICS_DBG("BackwardBranch in mol_branched " << std::endl);

	int b0 = first_b[generation];
	int bN = last_b[generation];
	std::vector<int> Br;
	std::vector<std::span<Real>> Gb;
	int M=lat->M;
	std::vector<Real> GS(4*M);
	auto gg_f = std::span<Real>(Gg_f);
	auto gg_b = std::span<Real>(Gg_b);
	auto unity_values = std::span<Real>(UNITY);
	int ss=0;
	for (int k = bN ; k >= b0 ; k--){
		if (Gnr[k]!=generation) {
			Br.clear(); Gb.clear();
			while (Gnr[k] != generation){
				Br.push_back(Gnr[k]);
				Gb.push_back(gg_f.subspan(static_cast<size_t>(last_s[Gnr[k]] * M), static_cast<size_t>(M)));
				ss=first_s[Gnr[k]];
				k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
			}
			Br.push_back(generation); ss--;
			Gb.push_back(gg_f.subspan(static_cast<size_t>(ss * M), static_cast<size_t>(M)));
			int length = Br.size();
			std::vector<Real> GX(length*M);
			for (int i=0; i<length; i++) std::copy_n(Gb[i].begin(), M, GX.data()+i*M);
			std::copy_n(gg_b.data()+((s+1)%2)*M, M, GS.data()+3*M);
			for (int i=0; i<length; i++) {
				std::copy_n(GS.data()+3*M, M, GS.data()+2*M);
				for (int j=0; j<length; j++) {
					if (i !=j) {
						std::copy_n(GX.data()+j*M, M, GS.begin());
						lat->propagate(GS.data(),unity_values.data(),0,1,M);
						for (int __i = 0; __i < M; ++__i) (GS.data()+2*M)[__i] = (GS.data()+2*M)[__i] * (GS.data()+M)[__i];
					}
				}
				std::copy_n(GS.data()+2*M, M, gg_b.begin());
				std::copy_n(GS.data()+2*M, M, gg_b.begin()+M);
				if (i<length-1) {
					BackwardBranch(Br[i],s);
				}
			}
			k++;
		} else {
			propagate_backward(Seg[mon_nr[k]]->G1.data(),s,k,M);
		}

	}
}

Real* mol_branched::ForwardBranch(int generation, int &s) {
NAMICS_DBG("ForwardBranch in mol_branched " << std::endl);
	int b0 = first_b[generation];
	int bN = last_b[generation];
	std::vector<int> Br;
	std::vector<std::span<Real>> Gb;
	int M=lat->M;
	std::vector<Real> GS(3*M);
	auto gg_f = std::span<Real>(Gg_f);
	auto unity_values = std::span<Real>(UNITY);

	Real* Glast=NULL;
	for (int k = b0; k<=bN ; ++k) {
		if (b0<k && k<bN) {
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1.data(),s,k,generation,M);
			} else {
				Br.clear(); Gb.clear();
				std::copy_n(Glast, M, GS.begin());
				while (Gnr[k] !=generation) {
					Br.push_back(Gnr[k]);
					Gb.push_back(std::span<Real>(ForwardBranch(Gnr[k],s), static_cast<size_t>(M)));
					k+=(last_b[Gnr[k]]-first_b[Gnr[k]]+1);
				}
				int length=Br.size();
				lat->propagate(GS.data(),Seg[mon_nr[k]]->G1.data(),0,2,M);

				for (int i=0; i<length; i++) {
					std::copy_n(Gb[i].begin(), M, GS.begin());
					lat->propagate(GS.data(),unity_values.data(),0,1,M);
					for (int __i = 0; __i < M; ++__i) (GS.data()+2*M)[__i] = (GS.data()+2*M)[__i] * (GS.data()+M)[__i];
				}
				std::copy_n(GS.data()+2*M, M, gg_f.begin()+s*M);
				s++;
			}
		} else {
			Glast=propagate_forward(Seg[mon_nr[k]]->G1.data(),s,k,generation,M);
		}
	}
	return Glast;
}



bool mol_branched::ComputePhi() {
NAMICS_DBG("ComputePhi in mol_branched " << std::endl);

	int M=lat->M;
	if (last_b.size() == 1) {
		// A single generation means an unbranched polymer, so the simple sweep is enough.
		int b0 = first_b[0];
		int bN = last_b[0];
		int s=0;
		Real* Glast=NULL;
		for (int b = b0; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1.data(),s,b,0,M);

		GN=lat->ComputeGN(Glast,M);

		s--;
		for (int b = bN ; b >= b0 ; b--) propagate_backward(Seg[mon_nr[b]]->G1.data(),s,b,M);
		return true;
	}

	int generation=0;
	int s=0;
	Real* G=ForwardBranch(generation,s);
	GN=lat->ComputeGN(G,M);
	s--;
	BackwardBranch(generation,s);

	return true;
}
