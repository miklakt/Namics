#include "molecule.h"
#include "mol_branched.h"


mol_branched::mol_branched(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_, std::string name_) : Molecule(In_,Lat_,Seg_,name_) {

}


mol_branched::~mol_branched() { }

void mol_branched::BackwardBra2ndO(Real* G_start, int generation,int &unity, int &s){
	(void)G_start;
NAMICS_DBG("BackwardBra2ndO in mol_branched " << std::endl);
	int b0 = first_b[generation];
	int bN = last_b[generation];
	std::vector<int> Br;
	std::vector<std::span<Real>> Gb;
	int M=lat->M;
	std::vector<Real> GS(3*M);
	std::vector<Real> GB(2*size*M);
	auto gg_f = std::span<Real>(Gg_f);
	auto gg_b = std::span<Real>(Gg_b);
	auto unity_values = std::span<Real>(UNITY);
	auto p = std::span<Real>(P);
	int ss=0;
	for (int k = bN ; k >= b0 ; k--){
		if (Gnr[k]!=generation) {
			Br.clear(); Gb.clear();
			while (Gnr[k] != generation){
				Br.push_back(Gnr[k]);
				Gb.push_back(gg_f.subspan(static_cast<size_t>(last_s[Gnr[k]] * M * size), static_cast<size_t>(M * size)));
				ss=first_s[Gnr[k]];
				k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
			}
			Br.push_back(generation); ss--;
			Gb.push_back(gg_f.subspan(static_cast<size_t>(ss * M * size), static_cast<size_t>(M * size))); //backbone last before branch point.
			int length = Br.size();
			std::vector<Real> GX(length*M);

			for (int i=0; i<length; i++) {
				lat->Terminate(GX.data()+i*M,Gb[i].data(),Markov,M);
			}

			std::copy_n(gg_b.data()+((s+1)%2)*M*size, M*size, GB.begin()); //Upto the branch point; no sides connected

			for (int i=0; i<length; i++) {
				std::copy_n(unity_values.begin(), M, GS.data()+2*M);
				for (int j=0; j<length; j++) {
					if (i !=j) {
						if (j==length-1) { //linking main chain
							std::copy_n(Gb[j].begin(), M*size, gg_b.begin());
							lat->propagateF(gg_b.data(),unity_values.data(),p.data(),0,1,M); //connect main chain including semiflexibility
							for (int __i = 0; __i < M * size; ++__i) GB[M*size + __i] = GB[__i] * gg_b[M*size + __i];
						} else { //linking sides
							std::copy_n(GX.data()+j*M, M, GS.begin());
							lat->propagate(GS.data(),unity_values.data(),0,1,M);
							for (int __i = 0; __i < M; ++__i) (GS.data()+2*M)[__i] = (GS.data()+2*M)[__i] * (GS.data()+M)[__i];
						}
					}
				}
				if (i<length-1) {
					for (int t=0; t<size; t++) for (int __i = 0; __i < M; ++__i) gg_b[t*M + __i] = GB[M*size + t*M + __i] * GS[2*M + __i]; //freely jointed onto branch
					std::copy_n(gg_b.begin(), M*size, gg_b.begin()+M*size);
					unity=-1;
					BackwardBra2ndO(gg_b.data(),Br[i],unity,s);

				} else { //prepare for main chain propagation
					for (int t=0; t<size; t++) for (int __i = 0; __i < M; ++__i) gg_b[t*M + __i] = GB[t*M + __i] * GS[2*M + __i]; //freely jointed onto main chain
					std::copy_n(gg_b.begin(), M*size, gg_b.begin()+M*size);
				}
			}
			k++;
		} else {
			propagate_backward(Seg[mon_nr[k]]->G1.data(),s,k,p.data(),unity,M);
			unity=0;
		}
	}
}


Real* mol_branched::ForwardBra2ndO(Real* G0, int generation, int &s) {
(void)G0;
NAMICS_DBG("ForwardBra2nd0 in mol_branched " << std::endl);
	int b0 = first_b[generation];
	int bN = last_b[generation];
	std::vector<int> Br;
	std::vector<std::span<Real>> Gb;
	int M=lat->M;
	std::vector<Real> GS(3*M);
	std::vector<Real> GB(2*size*M);
	auto gg_f = std::span<Real>(Gg_f);
	auto unity_values = std::span<Real>(UNITY);
	auto p = std::span<Real>(P);

	Real* Glast=NULL;
	for (int k = b0; k<=bN ; ++k) {
		if (b0<k && k<bN) {
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1.data(),s,k,p.data(),generation,M);
			} else {
				Br.clear(); Gb.clear();
				std::copy_n(Glast, M*size, GB.begin());
				while (Gnr[k] !=generation) { //collect information from branches.
					Br.push_back(Gnr[k]);
					Gb.push_back(std::span<Real>(ForwardBra2ndO(G0,Gnr[k],s), static_cast<size_t>(M * size)));
					k+=(last_b[Gnr[k]]-first_b[Gnr[k]]+1);
				}
				int length=Br.size();

				lat->propagateF(GB.data(),Seg[mon_nr[k]]->G1.data(),p.data(),0,1,M); //propagate main chain to branch point; keep semiflexibility

				std::copy_n(unity_values.begin(), M, GS.data()+2*M);
				for (int i=0; i<length; i++) {
					lat->Terminate(GS.data(),Gb[i].data(),Markov,M);
					lat->propagate(GS.data(),unity_values.data(),0,1,M);
					for (int __i = 0; __i < M; ++__i) (GS.data()+2*M)[__i] = (GS.data()+2*M)[__i] * (GS.data()+M)[__i];
				}
				for (int t=0; t<size; t++) for (int __i = 0; __i < M; ++__i) (GB.data()+M*size+t*M)[__i] = (GB.data()+M*size+t*M)[__i] * (GS.data()+2*M)[__i]; //all side freely jointed
				std::copy_n(GB.data()+M*size, M*size, gg_f.begin()+s*M*size);
				s++;
			}
		} else {
			Glast=propagate_forward(Seg[mon_nr[k]]->G1.data(),s,k,p.data(),generation,M);
		}
	}
	return Glast;
}



void mol_branched::BackwardBra(Real* G_start, int generation, int &s){
	(void)G_start;
NAMICS_DBG("BackwardBr in mol_branched " << std::endl);

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
					BackwardBra(G_start,Br[i],s);
				}
			}
			k++;
		} else {
			propagate_backward(Seg[mon_nr[k]]->G1.data(),s,k,generation,M);
		}

	}
}

Real* mol_branched::ForwardBra(Real* G0, int generation, int &s) {
(void)G0;
NAMICS_DBG("ForwardBra in mol_branched " << std::endl);
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
					Gb.push_back(std::span<Real>(ForwardBra(G0,Gnr[k],s), static_cast<size_t>(M)));
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
	bool success=true;
	int generation=0;
	int unity=0;
	int s=0;

	Real* G;
	if (Markov == 2) {
		G=ForwardBra2ndO(Seg[mon_nr[last_b[0]]]->G1.data(),generation,s);
	} else {
		G=ForwardBra(Seg[mon_nr[last_b[0]]]->G1.data(),generation,s);
	}
	GN=lat->ComputeGN(G,Markov,M);
	s--;
	if (Markov == 2) {
		BackwardBra2ndO(Seg[mon_nr[last_b[0]]]->G1.data(),generation,unity,s);
	} else {
		BackwardBra(Seg[mon_nr[last_b[0]]]->G1.data(),generation,s);
	}

	return success;
}
