#include "molecule.h"
#include "mol_branched.h"


mol_branched::mol_branched(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_, string name_) : Molecule(In_,Lat_,Seg_,name_) {

}


mol_branched::~mol_branched() { }

void mol_branched::BackwardBra2ndO(Real* G_start, int generation,int &unity, int &s){
	(void)G_start;
NAMICS_DBG("BackwardBra2ndO in mol_branched " << endl);
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=lat->M;
	Real* GS= (Real*) malloc(3*M*sizeof(Real));
	Real* GB= (Real*) malloc(2*size*M*sizeof(Real));
	int ss=0;
	for (int k = bN ; k >= b0 ; k--){
		if (Gnr[k]!=generation) {
			Br.clear(); Gb.clear();
			while (Gnr[k] != generation){
				Br.push_back(Gnr[k]);
				Gb.push_back(Gg_f+last_s[Gnr[k]]*M*size);
				ss=first_s[Gnr[k]];
				k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
			}
			Br.push_back(generation); ss--;
			Gb.push_back(Gg_f+ss*M*size); //backbone last before branch point.
			int length = Br.size();
			Real* GX= (Real*) malloc(length*M*sizeof(Real));

			for (int i=0; i<length; i++) {
				lat->Terminate(GX+i*M,Gb[i],Markov,M);
			}

			std::copy_n(Gg_b+((s+1)%2)*M*size, M*size, GB); //Upto the branch point; no sides connected

			for (int i=0; i<length; i++) {
				std::copy_n(UNITY, M, GS+2*M);
				for (int j=0; j<length; j++) {
					if (i !=j) {
						if (j==length-1) { //linking main chain
							std::copy_n(Gb[j], M*size, Gg_b);
							lat->propagateF(Gg_b,UNITY,P,0,1,M); //connect main chain including semiflexibility
							for (int __i = 0; __i < (M*size); ++__i) (GB+M*size)[__i] = (GB)[__i] * (Gg_b+M*size)[__i];
						} else { //linking sides
							std::copy_n(GX+j*M, M, GS);
							lat->propagate(GS,UNITY,0,1,M);
							for (int __i = 0; __i < (M); ++__i) (GS+2*M)[__i] = (GS+2*M)[__i] * (GS+M)[__i];
						}
					}
				}
				if (i<length-1) {
					for (int t=0; t<size; t++) for (int __i = 0; __i < (M); ++__i) (Gg_b + t*M)[__i] = (GB +M*size + t*M)[__i] * (GS+2*M)[__i]; //freely jointed onto branch
					std::copy_n(Gg_b, M*size, Gg_b+M*size);
					unity=-1;
					BackwardBra2ndO(Gg_b,Br[i],unity,s);

				} else { //prepare for main chain propagation
					for (int t=0; t<size; t++) for (int __i = 0; __i < (M); ++__i) (Gg_b + t*M)[__i] = (GB + t*M)[__i] * (GS+2*M)[__i]; //freely jointed onto main chain
					std::copy_n(Gg_b, M*size, Gg_b+M*size);
				}
			}
			free(GX);
			k++;
		} else {
			propagate_backward(Seg[mon_nr[k]]->G1,s,k,P,unity,M);
			unity=0;
		}
	}
	free(GS); free(GB);
}


Real* mol_branched::ForwardBra2ndO(Real* G0, int generation, int &s) {
NAMICS_DBG("ForwardBra2nd0 in mol_branched " << endl);
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=lat->M;
	Real* GS= (Real*) malloc(3*M*sizeof(Real));
	Real* GB= (Real*) malloc(2*size*M*sizeof(Real));

	Real* Glast=NULL;
	for (int k = b0; k<=bN ; ++k) {
		if (b0<k && k<bN) {
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,P,generation,M);
			} else {
				Br.clear(); Gb.clear();
				std::copy_n(Glast, M*size, GB);
				while (Gnr[k] !=generation) { //collect information from branches.
					Br.push_back(Gnr[k]);
					Gb.push_back(ForwardBra2ndO(G0,Gnr[k],s));
					k+=(last_b[Gnr[k]]-first_b[Gnr[k]]+1);
				}
				int length=Br.size();

				lat->propagateF(GB,Seg[mon_nr[k]]->G1,P,0,1,M); //propagate main chain to branch point; keep semiflexibility

				std::copy_n(UNITY, M, GS+2*M);
				for (int i=0; i<length; i++) {
					lat->Terminate(GS,Gb[i],Markov,M);
					lat->propagate(GS,UNITY,0,1,M);
					for (int __i = 0; __i < (M); ++__i) (GS+2*M)[__i] = (GS+2*M)[__i] * (GS+M)[__i];
				}
				for (int t=0; t<size; t++) for (int __i = 0; __i < (M); ++__i) (GB+M*size+t*M)[__i] = (GB+M*size+t*M)[__i] * (GS+2*M)[__i]; //all side freely jointed
				std::copy_n(GB+M*size, M*size, Gg_f+s*M*size);
				s++;
			}
		} else {
			Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,P,generation,M);
		}
	}
	free(GS); free(GB);
	return Glast;
}



void mol_branched::BackwardBra(Real* G_start, int generation, int &s){
	(void)G_start;
NAMICS_DBG("BackwardBr in mol_branched " << endl);

	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=lat->M;
	Real* GS= (Real*) malloc(4*M*sizeof(Real));
	int ss=0;
	for (int k = bN ; k >= b0 ; k--){
		if (Gnr[k]!=generation) {
			Br.clear(); Gb.clear();
			while (Gnr[k] != generation){
				Br.push_back(Gnr[k]);
				Gb.push_back(Gg_f+last_s[Gnr[k]]*M);
				ss=first_s[Gnr[k]];
				k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
			}
			Br.push_back(generation); ss--;
			Gb.push_back(Gg_f+ss*M);
			int length = Br.size();
			Real* GX= (Real*) malloc(length*M*sizeof(Real));
			for (int i=0; i<length; i++) std::copy_n(Gb[i], M, GX+i*M);
			std::copy_n(Gg_b+((s+1)%2)*M, M, GS+3*M);
			for (int i=0; i<length; i++) {
				std::copy_n(GS+3*M, M, GS+2*M);
				for (int j=0; j<length; j++) {
					if (i !=j) {
						std::copy_n(GX+j*M, M, GS);
						lat->propagate(GS,UNITY,0,1,M);
						for (int __i = 0; __i < (M); ++__i) (GS+2*M)[__i] = (GS+2*M)[__i] * (GS+M)[__i];
					}
				}
				std::copy_n(GS+2*M, M, Gg_b);
				std::copy_n(GS+2*M, M, Gg_b+M);
				if (i<length-1) {
					BackwardBra(G_start,Br[i],s);
				}
			}
			free(GX);
			k++;
		} else {
			propagate_backward(Seg[mon_nr[k]]->G1,s,k,generation,M);
		}

	}
	free(GS);
}

Real* mol_branched::ForwardBra(Real* G0, int generation, int &s) {
NAMICS_DBG("ForwardBra in mol_branched " << endl);
	int b0 = first_b[generation];
	int bN = last_b[generation];
	vector<int> Br;
	vector<Real*> Gb;
	int M=lat->M;
	Real* GS= (Real*) malloc(3*M*sizeof(Real));

	Real* Glast=NULL;
	for (int k = b0; k<=bN ; ++k) {
		if (b0<k && k<bN) {
			if (Gnr[k]==generation ){
				Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
			} else {
				Br.clear(); Gb.clear();
				std::copy_n(Glast, M, GS);
				while (Gnr[k] !=generation) {
					Br.push_back(Gnr[k]);
					Gb.push_back(ForwardBra(G0,Gnr[k],s));
					k+=(last_b[Gnr[k]]-first_b[Gnr[k]]+1);
				}
				int length=Br.size();
				lat->propagate(GS,Seg[mon_nr[k]]->G1,0,2,M);

				for (int i=0; i<length; i++) {
					std::copy_n(Gb[i], M, GS);
					lat->propagate(GS,UNITY,0,1,M);
					for (int __i = 0; __i < (M); ++__i) (GS+2*M)[__i] = (GS+2*M)[__i] * (GS+M)[__i];
				}
				std::copy_n(GS+2*M, M, Gg_f+s*M);
				s++;
			}
		} else {
			Glast=propagate_forward(Seg[mon_nr[k]]->G1,s,k,generation,M);
		}
	}
	free(GS);
	return Glast;
}



bool mol_branched::ComputePhi() {
NAMICS_DBG("ComputePhi in mol_branched " << endl);

	int M=lat->M;
	bool success=true;
	int generation=0;
	int unity=0;
	int s=0;

	Real* G;
	if (Markov == 2) {
		G=ForwardBra2ndO(Seg[mon_nr[last_b[0]]]->G1,generation,s);
	} else {
		G=ForwardBra(Seg[mon_nr[last_b[0]]]->G1,generation,s);
	}
	GN=lat->ComputeGN(G,Markov,M);
	s--;
	if (Markov == 2) {
		BackwardBra2ndO(Seg[mon_nr[last_b[0]]]->G1,generation,unity,s);
	} else {
		BackwardBra(Seg[mon_nr[last_b[0]]]->G1,generation,s);
	}

	return success;
}
