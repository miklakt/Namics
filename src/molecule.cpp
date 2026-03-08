#include "molecule.h"

Molecule::Molecule(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_, string name_) {
	In=In_; Seg=Seg_; name=name_;  Lat=Lat_;
NAMICS_DBG("Constructor for Mol " + name << endl);
	lat=Lat;
	KEYS.push_back("freedom");
	KEYS.push_back("composition");
	KEYS.push_back("ring");
	KEYS.push_back("theta");
	KEYS.push_back("phibulk");
	KEYS.push_back("n");
	KEYS.push_back("save_memory");
	KEYS.push_back("restricted_range");
	KEYS.push_back("compute_width_interface");
	KEYS.push_back("Markov");
	KEYS.push_back("k_stiff");
	KEYS.push_back("B");

	width=0;
	phi1=0;
	phiM=0;
	Dphi=0;
	pos_interface=0;
	ring=false;
	all_molecule=false;
	Markov =1;
	FillRangesList.clear();
	Filling=false;
	save_memory=false;
	J=0; Delta_MU=0; B=1;

}

Molecule::~Molecule() {
	DeAllocateMemory();
}

void Molecule :: DeAllocateMemory(){
NAMICS_DBG("DeallocateMemory for Mol " + name << endl);
	if (!all_molecule) return;
	free(H_phi);
	free(H_phitot);
	if (Markov==2) { free(P);}
	free(Gg_f);
	free(Gg_b);
	if (save_memory) free(Gs);
	free(UNITY);

	all_molecule=false;
}

void Molecule:: AllocateMemory() {
NAMICS_DBG("AllocateMemory in Mol " + name << endl);
	DeAllocateMemory();
	int M=lat->M;
	if (Markov==2){// && lat->lattice_type == simple_cubic) {
		int FJC = lat->FJC;
		P = (Real*) malloc(FJC*sizeof(Real));
		Real Q=0;
		KStiff=k_stiff;
		for (int k=0; k<FJC-1; k++) {
			P[k]=exp(-0.5*KStiff*(k*PIE/(FJC-1))*(k*PIE/(FJC-1)) ); //alternative to put U(theta)=-k(1-cos(theta))
			if (k>0) {
				if (lat->lattice_type==hexagonal) Q+= 2*P[k]; else Q+= 4*P[k]; //alternative is to use u_bend = Kstiff(1-cos(theta)), persistence length is l_p = b/ln <cos (theta)>
			} else Q=P[k];
		}
		P[FJC-1]=0; //Q+=P[FJC-1];
		if (lat->lattice_type==hexagonal&& !lat->stencil_full) Q*=2.0;
		for (int k=0; k<FJC; k++) { P[k]/=Q;
			cout << "P["<<k<<"] = " << P[k] << endl;
		}
	}

	if (save_memory) {
		int length_ = mon_nr.size();
		for (int i=0; i<length_; i++) last_stored.push_back(0);
		for (int i=0; i<length_; i++) {
			int n=int(pow(n_mon[i]*2,1.0/2.0)+0.5); //in sfbox apparently the pow 1/3 is reached. Here it fails for unknown reasons.
			if (n>n_mon[i]) n=n_mon[i]; //This seems to me to be enough...needs a check though..
			if (i==0) memory.push_back(n); else memory.push_back(n+memory[i-1]);
		}
	}
	N=0;
	if (save_memory) {
		N=memory[n_mon.size()-1];
	} else {
		int length_ = mon_nr.size();
		for (int i=0; i<length_; i++) {N+=n_mon[i];}
	}

	H_phi = (Real*) malloc(M*MolMonList.size()*sizeof(Real)); std::fill_n(H_phi, M * MolMonList.size(), 0);
	H_phitot = (Real*) malloc(M*sizeof(Real)); std::fill_n(H_phitot, M, 0);
	Gg_f = (Real*) malloc(M*N*sizeof(Real)*size);
	Gg_b = (Real*) malloc(M*2*sizeof(Real)*size);
	std::fill_n(Gg_f, M*N*size, 0);
	std::fill_n(Gg_b, 2*M*size, 0);
	phi=H_phi;
	rho=phi;
	if (save_memory) {Gs=(Real*) malloc(2*M*sizeof(Real)*size); std::fill_n(Gs, 2*M*size, 0);}
	phitot = H_phitot;
	UNITY = (Real*) malloc(M*sizeof(Real)*size); std::fill_n(UNITY, M*size, 0);
	all_molecule=true;
}

bool Molecule:: PrepareForCalculations(Real *KSAM) {
NAMICS_DBG("PrepareForCalculations in Mol " + name << endl);
int M=lat->M;
	std::copy_n(KSAM, M, UNITY);
	bool success=true;
	std::fill_n(phitot, M, 0);


		//lat->set_bounds(u+i*M);

		//lat->set_bounds(u+i*M);


		//lat->set_bounds(G1+i*M);
	//}
	std::fill_n(phi, M*MolMonList.size(), 0);


	return success;
}

bool Molecule::CheckInput(int start_, bool checking) {
NAMICS_DBG("Molecule:: CheckInput for mol " << name << endl);
start=start_;
phibulk=0;
n=0;
theta=0;
norm=0;
NAMICS_DBG("CheckInput for Mol " + name << endl);
	bool success=true;
	if (!In->CheckParameters("mol",name,start, KEYS, PARAMETERS)) {
		success=false;
	} else {
		save_memory=false;
		if (GetValue("save_memory").size()>0) {
			save_memory=ParseBool(GetValue("save_memory"),false);
		}
		if (GetValue("composition").size()==0) {cout << "For mol '" + name + "' the definition of 'composition' is required" << endl; success = false;
		} else {
			try {
				if (!Decomposition(GetValue("composition"))) {
					cout << "For mol '" + name + "' the composition is rejected. " << endl;
					success=false;
				}
			} catch (const char* error) {
				cerr << error << endl;
				success = false;
			}
		}
		if (GetValue("restricted_range").size()>0) {
			if (GetValue("freedom")!="range_restricted") cout <<"For mol '" + name + "' freedom is not set to 'range_restricted' and therefore  the value of 'restricted_range' is ignored" << endl;
		}

		if (checking) {
			if (GetValue("freedom").size() > 0) freedom = GetValue("freedom"); else {
				cout <<"For molecule " << name << " no value for 'freedom' was found " << endl;
				success=false;
			}
			return success; //because moltype and freedom are known; as start <0 the checkinput can be terminated.
		}
		if (IsPinned()) {
			if (GetValue("freedom").size()==0) {
					cout <<"For mol " + name + " the setting for 'freedom' was not set" << endl; return false;
			} else {
				vector<string> free_list;
				free_list.push_back("restricted");
				free_list.push_back("fill_range");
				if (!ParseString(GetValue("freedom"),freedom,free_list,"In mol " + name + " the value for 'freedom' is not recognised ")) return false;
				if (freedom=="restricted") {
					if (GetValue("theta").size() ==0 && GetValue("n").size()==0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' or 'freedom = range_restricted',should be combined with a value for 'theta' or 'n'; do not use both settings! "<<endl; success=false;
					} else {
							if (GetValue("theta").size() >0 && GetValue("n").size()>0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' of 'freedom = range_restricted' do not specify both 'n' and 'theta' "<<endl; success=false;
					} else {
							if (GetValue("n").size()>0) {n=ParseReal(GetValue("n"),10*lat->volume);theta=n*chainlength;}
							if (GetValue("theta").size()>0) {theta = ParseReal(GetValue("theta"),10*lat->volume);n=theta/chainlength;}
							if (theta < 0 ) {    //|| theta > lat->volume) {
								cout << "In mol " + name + ", the value of 'n' or 'theta' " << theta << "  is out of range 0 .. 'volume'/N, cq 'volume' "<< lat->volume << endl; success=false;
							}
						}
					}
				}
				if (freedom=="fill_range") {
					freedom="restricted";
					Filling=true;
					if (GetValue("theta").size() > 0 || GetValue("n").size()>0 || GetValue("phibulk").size() > 0) {
						if (start==1) cout <<"For mol " + name + " the freedom is set to 'fill-range-of{-mon_name}' and therefore the value of 'theta', the value of 'n', or the value of 'phibulk' is ignored. " << endl;
					}
				}
			}
		} else
		if (GetValue("freedom").size()==0 && !IsTagged() ) {
			cout <<"For mol " + name + " the setting 'freedom' is expected: options: 'free' 'restricted' 'solvent' 'neutralizer' 'range_restricted' 'tagged' . Problem terminated " << endl; success = false;
			} else {

				if (!IsTagged()) {
				vector<string> free_list;
				if (!IsPinned()) {
					free_list.push_back("free");
					free_list.push_back("solvent");
					free_list.push_back("neutralizer");
					free_list.push_back("range_restricted");
				}
				free_list.push_back("restricted");
				if (!ParseString(GetValue("freedom"),freedom,free_list,"In mol " + name + " the value for 'freedom' is not recognised ")) success=false;
				if (freedom == "solvent") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the solvent" << endl; }
				}
				if (freedom == "neutralizer") {
					if (IsPinned()) {success=false; cout << "Mol '" + name + "' is 'pinned' and therefore this molecule can not be the neutralizer" << endl; }
					if (!IsCharged()) {success=false; cout << "Mol '" + name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << endl; }
				}
				if (freedom == "free") {
					if (GetValue("phibulk").size() ==0) {
						cout <<"In mol " + name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<endl; return false;
					} else {
						phibulk=ParseReal(GetValue("phibulk"),-1);
						if (phibulk < 0 || phibulk >1) {
							cout << "In mol " + name + ", the value of 'phibulk' is out of range 0 .. 1." << endl; return false;
						}
					}
				}

				B=1;
				if (GetValue("B").size()>0){
					B=ParseReal(GetValue("B"),B);
					if (B<1e-9) {
						cout <<"for Mol" + name + " mobility B should have a posititve value. Default value B=1 is chosen. " << endl;
						B=1;
					}
				}

				if (freedom == "restricted" || freedom=="range_restricted") {
					//} else {
						if (GetValue("theta").size() ==0 && GetValue("n").size()==0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' or 'freedom = range_restricted',should be combined with a value for 'theta' or 'n'; do not use both settings! "<<endl; success=false;
						} else {
							if (GetValue("theta").size() >0 && GetValue("n").size()>0) {
							cout <<"In mol " + name + ", the setting 'freedom = restricted' of 'freedom = range_restricted' do not specify both 'n' and 'theta' "<<endl; success=false;
							} else {

								if (GetValue("n").size()>0) {n=ParseReal(GetValue("n"),10*lat->volume);theta=n*chainlength;}
								if (GetValue("theta").size()>0) {theta = ParseReal(GetValue("theta"),10*lat->volume);n=theta/chainlength;}
								if (theta < 0 || theta > lat->volume) {
									cout << "In mol " + name + ", the value of 'n' or 'theta' is out of range 0 .. 'volume', cq 'volume'/N." << endl; success=false;

								}
							}
						}
					//}
				}
				if (freedom =="range_restricted" ) {
					if (GetValue("restricted_range").size() ==0) {
						success=false;
						cout<<"In mol '" + name + "', freedom is set to 'range_restricted'. In this case we expect the setting for 'restricted_range'. This setting was not found. Problem terminated. " << endl;
					} else { //read range;
						int *HP=NULL;
						int M=lat->M;
						int npos=0;
						bool block;
						R_mask=(Real*)malloc(M*sizeof(Real));
						string s="restricted_range";
						int *r=(int*) malloc(6*sizeof(int));
						success=lat->ReadRange(r,HP,npos,block,GetValue("restricted_range"),0,name,s);
						lat->CreateMASK(R_mask,r,HP,npos,block);
						theta_range = theta;
						n_range = theta_range/chainlength;
						free(r);
					}

				}

			} else {
				if (GetValue("theta").size() >0 || GetValue("n").size() > 0 || GetValue("phibulk").size() >0 || GetValue("freedom").size() > 0) cout <<"Warning. In mol " + name + " tagged segment(s) were detected. In this case no value for 'freedom' is needed, and also 'theta', 'n' and 'phibulk' values are ignored. " << endl;
			}
		}


		ring=false;
		if (GetValue("ring").size() > 0) {
				ParseBool(GetValue("ring"),ring,"Input for ring is either 'true' or 'false'. Moreover, first and last segments of the backbone will be put on top of each other (chain length gets shorter by one). ");
				if (ring) {
					int length;
					switch (MolType) {
					case monomer:
						cout << "Can not make a ring from a molecule type monomer; ring option is ignored" << endl;
						ring=false;
						break;
					case linear:
					case branched:
						length=mon_nr.size();
						if (mon_nr[0]==mon_nr[length-1]) {
							if (n_mon[0] !=1 || n_mon[length-1] !=1) {
								cout <<"Fist and last block of main chain should consist of just one segment: e.g. (A)1(B)99(A)1 is a ring of 100 segments (one A only)."  << endl;
								ring =false;
							} else  {
								cout <<"ring can be implemented " << endl;
								cout <<"chain length reduced by one" << endl;
								chainlength--;
							}
						} else {
							ring =false;
							cout <<"first and last segment of the main chain should be of the same type, as these two are 'merged'. Chain length is reduced by one." << endl;
							cout <<"make sure to start the main chain and end the main chain by a block with length 1: for example '(B)1(B)99(B)1' will be a ring of 100 segments." << endl;
						}
					break;
					default:
						ring=false;
					break;

				}
			}
		}

	}
	Markov=1;
	if (GetValue("Markov").size()>0) Markov=ParseInt(GetValue("Markov"),1);
	if (Markov<1 || Markov>2) {
		cout <<" Integer value for 'Markov' is by default 1 and may be set to 2 for some mol_types and fjc-choices only. Markov value out of bounds. Proceed with caution. " << endl; success = false;
	}
	if (Markov==2) lat->Markov=2;
	k_stiff=lat->k_stiff; //pick up 'default' value from lattice.
	if (GetValue("k_stiff").size()>0) {
		k_stiff=ParseReal(GetValue("k_stiff"),k_stiff);
		if (k_stiff<0 || k_stiff>10) {
			success =false;
			cout <<" Real value for 'k_stiff' out of bounds (0 < k_stiff < 10). " << endl;
			cout <<" For Markov == 2: u_bend (theta) = 0.5 k_stiff theta^2, where 'theta' is angle for bond direction deviating from the straight direction. " <<endl;
			cout <<" You may interpret 'k_stiff' as the molecular 'persistence length' " << endl;
			cout <<" k_stiff is a 'default value'. Use molecular specific values to overrule the default when appropriate (future implementation....) " << endl;
		}
		if (lat->fjc>1 && lat->gradients>1) {
			success=false;
			cout <<" Work in progress.... Currently, Markov == 2 is implemented in gradients>1 for FJC_choices = 3 " << endl;
		}
	}
	size=0;
	if (Markov ==2) {
		if (lat->gradients==1) size = lat->FJC;
		if (lat->gradients==2) { //assume fjc=1...
		    if (lat->stencil_full) {
				lat->stencil_full=false; cout <<"Warning: stencil_full is set to false" << endl;
			}
			if (lat->lattice_type==hexagonal ) size =  12;
			if (lat->lattice_type==simple_cubic ) size = 2*lat->FJC-1;
			if (lat->lattice_type==hexagonal) { success=false;
				cout <<"Warning: Markov = 2 & gradients=2 lattice_type = hexagonal...not certified. Caution recommended, even without obvious error messages..." << endl;
			}
		}
		if (lat->gradients==3) {
		    if (lat->stencil_full) {
				lat->stencil_full=false; cout <<"Warning: stencil_full is set to false" << endl;
			}

			if (lat->lattice_type==hexagonal) size=12;
			if (lat->lattice_type==simple_cubic) size = 6;
			if (lat->lattice_type==hexagonal) { success=false;
				cout <<"Warning: Markov = 2 & gradients=3 lattice_type = hexagonal...not certified!!!. Caution recommended, even without obvious error messages...." << endl;
			}
		}
	} else size = 1;
	if (size==0) {success=false; cout <<"Attention: size in molecule is not set; combination gradients (1,2,3),  Markov=2, stencil_full (true,false), lattice_type (hexagonal, simple_cubic) not implemented" << endl;}

	return success;
}

int Molecule::GetMonNr(string s){
NAMICS_DBG("GetMonNr for Mon " + name << endl);
	int n_segments=In->MonList.size();
	int found=-1;
	int i=0;
	while(i<n_segments) {
		if (Seg[i]->name ==s) found=i;
		i++;
	}
	return found;
}

bool Molecule::ExpandBrackets(string &s) {
NAMICS_DBG("Molecule:: ExpandBrackets" << endl);
	bool success=true;
	if (s[0] != '(') {cout <<"illegal composition. Expects composition to start with a '(' in: " << s << endl; return false;}
	vector<int> open;
	vector<int> close;
	bool done=false; //now interpreted the (expanded) composition
	while (!done) { done = true;
		open.clear(); close.clear();
		if (!In->EvenBrackets(s,open,close)) {
			cout << "s : " << s << endl;
			cout << "In composition of mol '" + name + "' the backets are not balanced."<<endl; success=false; return success;
		 }
		int length=open.size();
		int pos_open;
		int pos_close;
		int pos_low=0;
		int i_open=0; {pos_open=open[0]; pos_low=open[0];}
		int i_close=0; pos_close=close[0];
		if (pos_open > pos_close) {cout << "Brackets open in composition not correct" << endl; return false;}
		while (i_open < length-1 && done) {
			i_open++;
			pos_open=open[i_open];
			if (pos_open < pos_close && done) {
				i_close++; if (i_close<length) pos_close=close[i_close];
			} else {
				if (pos_low==open[i_open-1] ) {
					pos_low=pos_open;
					i_close++; if (i_close<length) pos_close=close[i_close];
					if (pos_open > pos_close) {cout << "Brackets open in composition not correct" << endl; return false;}
				} else {
					done=false;

					int x=ParseInt(s.substr(pos_close+1),-1);
						if (x<1) {
								cout <<"Number of repeats must be a positive integer in composition at pos : " << pos_close+1 << " for: " << s << endl; return false;
						}
					string sA,sB,sC;
					if (s.substr(pos_open-1,1)=="]") {pos_open --;  }
					sA=s.substr(0,pos_low);
					sB=s.substr(pos_low+1,pos_close-pos_low-1);
					sC=s.substr(pos_open,s.size()-pos_open+1);
					s=sA;for (int k=0; k<x; k++) s.append(sB); s.append(sC);
				}
			}
		}
		if (pos_low < open[length-1]&& done) {
			done=false;
			pos_close=close[length-1];
			int x=ParseInt(s.substr(pos_close+1),0);
			string sA,sB,sC;
			sA=s.substr(0,pos_low);
			sB=s.substr(pos_low+1,pos_close-pos_low-1);
			sC="";
			s=sA;for (int k=0; k<x; k++) s.append(sB); s.append(sC);
		}

	}
	if (s[s.size()-1]==']') {cout <<"illegal composition. Composition can not end with a ']' in: " << s << endl; return false;}
	return success;
}

bool Molecule::Interpret(string s,int generation){
NAMICS_DBG("Molecule:: Interpret" << endl);
	if (s=="[") return true;
	bool success=true;
	vector<int>open;
	vector<int>close;
	In->EvenBrackets(s,open,close);
	if (open.empty()) {
		cout << "In composition of mol '" + name + "' an invalid token was found: " << s << endl;
		return false;
	}
	int k=0;
	int length=open.size();
	while (k<length) {
		string segname=s.substr(open[k]+1,close[k]-open[k]-1);
		int mnr=GetMonNr(segname);
		if (mnr <0)  {cerr <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised"  << endl; success = false;
		throw "Composition Error";
		} else {

			int stored=Gnr.size();
			if (stored>0) {//fragments at branchpoint need to be just 1 segment long.
				if (Gnr[stored-1]<generation) {
					if (n_mon[stored-1]>1) {
						n_mon[stored-1]--;
						n_mon.push_back(1);
						mon_nr.push_back(mon_nr[stored-1]);
						Gnr.push_back(Gnr[stored-1]);
					 	last_b[Gnr[stored-1]]++;
					}
				}
			}
			mon_nr.push_back(mnr);
			Gnr.push_back(generation);
			if (first_s[generation] < 0) first_s[generation]=chainlength;
			if (first_b[generation] < 0) first_b[generation]=mon_nr.size()-1;
			last_b[generation]=mon_nr.size()-1;
		}
		int nn = ParseInt(s.substr(close[k]+1,s.size()-close[k]-1),0);
		if (nn<1) {cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity " << endl; success=false; return success;
		} else {
			n_mon.push_back(nn);
		}
		chainlength +=nn; last_s[generation]=chainlength;
		k++;
	}
	return success;
}

bool Molecule::GenerateTree(string s,int generation,int &pos, vector<int> open,vector<int> close) {
NAMICS_DBG("Molecule:: GenerateTree" << endl);
	bool success=true;
	string ss;
	int i=0;
	int newgeneration=0;
	int new_generation=0;
	int length=open.size();
	int pos_open=0;
	int pos_close=s.length();
	bool openfound,closedfound;
	while  (pos_open<pos_close && success) {
		pos_open =s.length()+1;
		pos_close=s.length();
		openfound=closedfound=false;
		i=0;
		while (i<length && !(openfound && closedfound) ){

			if (close[i]>pos && !closedfound) {closedfound=true; pos_close=close[i]+1; new_generation=i+1;}
			if (open[i]>=pos && !openfound) {openfound=true; pos_open=open[i]+1; newgeneration=i+1;}
			i++;
		}

		if (pos_close<pos_open) {
			ss=s.substr(pos,pos_close-pos);
			if (ss.substr(0,1)=="[") {
				pos=pos+1;
				first_s.push_back(-1);
				last_s.push_back(-1);
				first_b.push_back(-1);
				last_b.push_back(-1);
				success=GenerateTree(s,new_generation,pos,open,close);
				if (!success) {cout <<"error in generate tree ." << endl; return success; }
				pos_close=pos_open+1;
			} else {
				pos=pos_close;
				success =Interpret(ss,generation);
				if (!success)  {cout <<"error in interpret ." << endl; return success; }
			}
		} else {
			ss=s.substr(pos,pos_open-pos);
			pos=pos_open;
			success =Interpret(ss,generation);
			if (!success) {cout <<"error in Interpret" << endl;  return success;}
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			success=GenerateTree(s,newgeneration,pos,open,close);
			if (!success) {cout <<"error in generate tree .." << endl; return success;}
		}
	}
	return success;
}

//}

bool Molecule::Decomposition(string s){
NAMICS_DBG("Decomposition for Mol " + name << endl);
	bool success = true;
	MolType=linear;
	chainlength=0;
	vector<int> open;
	vector<int> close;

	if (!ExpandBrackets(s)) {success=false; return success;}
	if (!In->EvenSquareBrackets(s,open,close)) {
		cout << "Error in composition of mol '" + name + "'; the square brackets are not balanced in " << s << endl;
		success=false; return success;
	}
	if (open.size()>0) {
		MolType=branched;
	}
	int generation=0;
	int pos=0;
	MolMonList.clear();
	int lopen;
	first_s.clear();
	last_s.clear();
	first_b.clear();
	last_b.clear();
	first_s.push_back(-1);
	last_s.push_back(-1);
	first_b.push_back(-1);
	last_b.push_back(-1);
	if (MolType==branched) {
		lopen=open.size();
		for (int i=1; i<lopen; i++) {
			if ((open[i]-open[i-1])==1 || (close[i]-close[i-1])==1) {
				cout <<"In molecule " + name + " in 'composition', two similar square brackets in a row '[[' or ']]' is not allowed" << endl;
				success=false;
				return success;
			}
		}
	}
	success = GenerateTree(s,generation,pos,open,close);
	if (!success) {
		cout << "GenerateTree failed" << endl;
		return success;
	}
	if (MolType==branched) { //invert numbers;
		int g_length=first_s.size();
		int length = n_mon.size();
		int xxx;
		int ChainLength=last_s[0];

		for (int i=0; i<g_length; i++) {
			first_s[i] = ChainLength-first_s[i]-1;
			last_s[i] = ChainLength-last_s[i]-1;
			xxx=first_s[i]; first_s[i]=last_s[i]+1; last_s[i]=xxx;
			first_b[i]=length-first_b[i]-1;
			last_b[i]=length-last_b[i]-1;
			xxx=first_b[i]; first_b[i]=last_b[i]; last_b[i]=xxx;
		}

		for (int i=0; i<length/2; i++) {
			xxx=Gnr[i]; Gnr[i]=Gnr[length-1-i]; Gnr[length-1-i]=xxx;
			xxx=n_mon[i]; n_mon[i]=n_mon[length-1-i]; n_mon[length-1-i]=xxx;
			xxx=mon_nr[i]; mon_nr[i]=mon_nr[length-1-i]; mon_nr[length-1-i]=xxx;
		}
		//		}
		//}

	}

	success=MakeMonList();
	if (chainlength==1) MolType=monomer;
	return success;
}

int Molecule::GetChainlength(void){
NAMICS_DBG("GetChainlength for Mol " + name << endl);
	return chainlength;
}

bool Molecule:: MakeMonList(void) {
NAMICS_DBG("Molecule:: MakeMonList" << endl);
	MolMonList.clear();
	bool success=true;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (!In->InSet(MolMonList,mon_nr[i])) {
			if (Seg[mon_nr[i]]->GetFreedom()=="frozen") {
				success = false;
				cout << "In 'composition of mol " + name + ", a segment was found with freedom 'frozen'. This is not permitted. " << endl;

			}
			MolMonList.push_back(mon_nr[i]);
		}
		i++;
	}
	i=0;
	int pos;
	while (i<length) {
		if (In->InSet(MolMonList,pos,mon_nr[i])) {molmon_nr.push_back(pos);
		} else {cout <<"program error in mol PrepareForCalcualations" << endl; }
		i++;
	}
	return success;
}

bool Molecule::IsPinned() {
NAMICS_DBG("IsPinned for Mol " + name << endl);
	bool success=false;
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
        	if (Seg[MolMonList[i]]->GetFreedom()=="pinned") success=true;
		i++;
	}
	return success;
}

int Molecule::GetPinnedSeg() {
NAMICS_DBG("GetPinnedSeg for Mol " + name << endl);
	int segnr=-1;
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
        	if (Seg[MolMonList[i]]->GetFreedom()=="pinned") {
				if (segnr > -1) cout <<"multiple segnrs in GetPinnedSeg() " << endl; else segnr=MolMonList[i];
			}
		i++;
	}
	return segnr;
}

bool Molecule::IsTagged() {
NAMICS_DBG("IsTagged for Mol " + name << endl);
	bool success=false;
	int length=MolMonList.size();
	int i=0;
	while (i<length) {
		if (Seg[MolMonList[i]]->freedom=="tagged") {success = true; tag_segment=MolMonList[i]; }
		i++;
	}
	return success;
}

Real Molecule::Charge() {
NAMICS_DBG("Molecule:: Charge" << endl);
	Real charge=0;
	int length=mon_nr.size();
	int length_states;
	for (int i=0; i<length; i++) {
		if (Seg[mon_nr[i]]->state_name.size() >1) {
			length_states=Seg[mon_nr[i]]->state_name.size();
			for (int j=0; j<length_states; j++) charge +=Seg[mon_nr[i]]->state_alphabulk[j]*Seg[mon_nr[i]]->state_valence[j]*n_mon[i];
		} else	charge +=Seg[mon_nr[i]]->valence*n_mon[i];
	}
	return charge/chainlength;
}

bool Molecule::IsCharged() {
NAMICS_DBG("IsCharged for Mol " + name << endl);
	Real charge =0;
	bool ischarged=false;
	int length = n_mon.size();
	int length_states;
	int i=0;
	while (i<length) {
		if (Seg[mon_nr[i]]->state_name.size()>0) {
			length_states=Seg[mon_nr[i]]->state_name.size();
			for (int j=0; j<length_states; j++) {if (Seg[mon_nr[i]]->state_valence[j]!=0) ischarged=true; }
		} else
		charge +=n_mon[i]*Seg[mon_nr[i]]->valence;
		i++;
	}
	if (charge !=0) ischarged=true;
	return ischarged;
}

void Molecule::PutParameter(string new_param) {
NAMICS_DBG("PutParameter for Mol " + name << endl);
	KEYS.push_back(new_param);
}

string Molecule::GetValue(string parameter) {
	auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

void Molecule::push(string s, Real X) {
NAMICS_DBG("push (Real) for Mol " + name << endl);
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Molecule::push(string s, int X) {
NAMICS_DBG("push (int) for Mol " + name << endl);
	ints.push_back(s);
	ints_value.push_back(X);
}
void Molecule::push(string s, bool X) {
NAMICS_DBG("push (bool) for Mol " + name << endl);
	bools.push_back(s);
	bools_value.push_back(X);
}
void Molecule::push(string s, string X) {
NAMICS_DBG("push (string) for Mol " + name << endl);
	strings.push_back(s);
	strings_value.push_back(X);
}



void Molecule::PushOutput() {
NAMICS_DBG("PushOutput for Mol " + name << endl);
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("composition",GetValue("composition"));
	if (IsTagged()) {string s="tagged"; push("freedom",s);} else {push("freedom",freedom);}
	if (freedom=="free") theta = lat->WeightedSum(phitot);
	push("Markov",Markov);
	push("k_stiff",k_stiff);
	if (lat->gradients==3) {
		int MZ=lat->MZ;
		int MY=lat->MY;
		int MX=lat->MX;
		int JX=lat->JX;
		int JY=lat->JY;
		for (int z=1; z<MZ+1; z++) {
			Real phiz=0;
			for (int x=1; x<MX+1; x++) for (int y=1;y<MY+1;y++) {
				phiz +=H_phitot[x*JX+y*JY+z];
			}
			phiz /= MX*MY;
			if (z==1) push("phiz[1]",phiz);
			if (z==2) push("phiz[2]",phiz);
			if (z==3) push("phiz[3]",phiz);
			if (z==4) push("phiz[4]",phiz);
			if (z==5) push("phiz[5]",phiz);
			if (z==6) push("phiz[6]",phiz);
			if (z==7) push("phiz[7]",phiz);
			if (z==8) push("phiz[8]",phiz);
			if (z==9) push("phiz[9]",phiz);
			if (z==10) push("phiz[10]",phiz);
			if (z==11) push("phiz[11]",phiz);
			if (z==12) push("phiz[12]",phiz);
			if (z==13) push("phiz[13]",phiz);
			if (z==14) push("phiz[14]",phiz);
			if (z==15) push("phiz[15]",phiz);
			if (z==16) push("phiz[16]",phiz);
			if (z==17) push("phiz[17]",phiz);
			if (z==18) push("phiz[18]",phiz);
			if (z==19) push("phiz[19]",phiz);
			if (z==20) push("phiz[20]",phiz);
		}
	}
	if (Markov==2) {
		for (int k=0; k<size; k++){
			if (k==0) push("P[0]",P[0]);
			if (k==1) push("P[1]",P[1]);
			if (k==2) push("P[2]",P[2]);
			if (k==3) push("P[3]",P[3]);
			if (k==4) push("P[4]",P[4]);
		}
	}
	//	    lat->remove_bounds(phitot);
	//	    (theta) = 0; for (int __i = 0; __i < (lat->M); ++__i) (theta) += (phitot)[__i];
	//  }
	//}

	push("Rg",pow((lat->Moment(phitot,0.0,2)/chainlength),0.5));
	lat->remove_bounds(phitot);
	theta=lat->WeightedSum(phitot);
	push("theta",theta);
	Real thetaexc=theta-lat->volume*phibulk;
	push("theta_exc",thetaexc);
	push("n_exc",thetaexc/chainlength);
	push("nexc",thetaexc/chainlength);
	push("thetaexc",thetaexc);
	push("theta_Gibbs",theta_Gibbs);
	if (R_Gibbs>0) push("R_Gibbs",R_Gibbs);
	push("n",n);
	push("chainlength",chainlength);
	push("phibulk",phibulk);
	push("Mu",Mu);
	push("mu",Mu); push("MU",Mu);
	if (lat->gradients==3) {
		Real TrueVolume=lat->MX*lat->MY*lat->MZ;
		Real Volume_particles=0;
		int num_of_seg=In->MonList.size();
		for (int i=0; i<num_of_seg; i++) {
			if (Seg[i]->freedom=="frozen") {
				for (int __j = 0; __j < lat->M; ++__j) Volume_particles += Seg[i]->MASK[__j];
			}
		}
		push("Gamma",theta-(TrueVolume-Volume_particles)*phibulk);
	}
	if (GetValue("compute_width_interface").size()>0){
		if (!ComputeWidth()) {
			cout <<"Computation of width of interface is rejected" <<endl;
		}
	}
	push("width",width);
	push("phi1",H_phitot[lat->fjc]);
	push("phiM",H_phitot[lat->M-2*lat->fjc]);
	push("Dphi",phi1-phiM);
	push("pos_interface",pos_interface);
	push("phi_average",phi_av);
	if (chainlength==1) {
		int seg=MolMonList[0];
		if (Seg[seg]->ns >1) {
			if (mu_state.size() ==0) for (int i=0; i<Seg[seg]->ns; i++) mu_state.push_back(Mu);
			for (int i=0; i<Seg[seg]->ns; i++) {
				mu_state[i]+=log(Seg[seg]->state_alphabulk[i]);
				push("mu-"+Seg[seg]->state_name[i],mu_state[i]);
			}
		}
	}
	int M=lat->M;
	Real phimax=H_phitot[M/2];
	bool maxfound=false;
	int i=M/2;
	while (!maxfound) {
		i++;
		if (H_phitot[i]> phimax ) phimax =H_phitot[i]; else maxfound=true;
	}
	maxfound=false;
	i=M/2;
	while (!maxfound) {
		i--;
		if (H_phitot[i]> phimax ) phimax =H_phitot[i]; else maxfound=true;
	}

	push("phiMax",phimax);

	push("GN",GN);
	push("norm",norm);
	J=0;

	int molmonlength=MolMonList.size();
	for (int i=0; i<molmonlength; i++) {
		J+=Seg[MolMonList[i]]->J;
	}
	push("J",J/chainlength);
	push("DeltaMu",Delta_MU);
	string s="profile;0"; push("phi",s);
	int length = MolMonList.size();
	for (int i=0; i<length; i++) {
		stringstream ss; ss<<i+1; string str=ss.str();
		s= "profile;"+str; push("phi_"+Seg[MolMonList[i]]->name,s);
	}
}

Real* Molecule::GetPointer(string s, int &SIZE) {
NAMICS_DBG("GetPointer for Mol " + name << endl);
	vector<string> sub;
	int M= lat->M;
	In->split(s,';',sub);
	if (sub[0]=="profile") {
		SIZE=M;

		if (sub[1]=="0") {
			lat->set_bounds(phitot);
			return H_phitot;
		}

		int length=MolMonList.size();
		int i=0;
		while (i<length) {
			stringstream ss; ss<<i+1; string str=ss.str();
			if (sub[1]==str) {
				lat->set_bounds(phi+i*M);
				return H_phi+i*M;
			}
			i++;
		}
	}
	return NULL;
}

void Molecule::PutTheta(Real T){
	theta=T;
	n=theta/chainlength;
}

int* Molecule::GetPointerInt(string s, int &SIZE) {
	(void)SIZE;
NAMICS_DBG("GetPointerInt for Mol " + name << endl);
	vector<string> sub;
	In->split(s,';',sub);
	if (sub[0]=="array") { //set SIZE and return array pointer.
	}
	return NULL;
}
int Molecule::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
NAMICS_DBG("GetValue (long) for Mol " + name << endl);
	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = Reals.size();
	while (i<length) {
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true"; else string_result="false";
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}
		i++;
	}
	return 0;
}

Real Molecule::ComputeGibbs(Real R_gibbs) {
NAMICS_DBG("ComputeGibbs for Mol " + name << endl);
	int fjc=lat->fjc;
	int M=lat->M;
	int gradients=lat->gradients;
	if (gradients>1) {cout <<"Error in ComputeGibbs; gadients not equal to 1 " << endl; return 0.0;}
	lat->remove_bounds(phitot);
	Real phi_low =phitot[fjc+1];
	Real phi_high=phitot[M-2*fjc-1];
	theta=lat->WeightedSum(phitot);
	Real theta_exc=theta-lat->volume*phibulk;
	if (freedom=="solvent") {
		theta_Gibbs=0;
		R_Gibbs=theta_exc/(phi_low-phi_high);
		return  R_Gibbs;
	} else {
		R_Gibbs=R_gibbs;
		theta_Gibbs=theta_exc-(phi_low-phi_high)*R_gibbs;
		return theta_Gibbs;
	}
}

bool Molecule::ComputeWidth() {
NAMICS_DBG("ComputeWidth for Mol " + name << endl);
	bool success=true;
	int M=lat->M;
	if (lat->gradients>1) {success=false; cout <<" Compute width of interface only in system with 'one-gradient'" << endl; return success; }
	if (GetValue("compute_width_interface")!="true") {
		cout <<"Interfacial width not computed because value for 'compute_width_interface' was not set to 'true'. " << endl; success=false; return success;
	} else {
		int fjc=lat->fjc;
		width=0;
		Dphi=0;
		phi_av=0;
		phi1=phitot[fjc]; phiM=phitot[M-2*fjc]; Dphi=phi1-phiM;
		if (Dphi*Dphi<1e-20) {cout << "There is no interface present. Width evaluation failed;" << endl; return success; }
		for (int x=fjc; x<M-fjc-1; x++) {
			if ((phitot[x]-phitot[x+1])/Dphi > width) {width = (phitot[x]-phitot[x+1])/Dphi; pos_interface=x+0.5; phi_av=(phitot[x]+phitot[x+1])/2;}
			}
	}
	if (width >0) width = 1.0/width/lat->fjc;
	pos_interface = (pos_interface)/lat->fjc;
	return success;
}

void Molecule::NormPerBlock(int split) {
	int MX=lat->MX/split;
	int MY=lat->MY/split;
	int MZ=lat->MZ/split;
	int JX=lat->JX;
	int JY=lat->JY;
	int M = lat->M;
	Real theta_block;
	int blocknr=-1;
	for(int i=0; i<split; i++)
	for(int j=0; j<split; j++)
	for(int k=0; k<split; k++){
		theta_block=0;
		blocknr++;
		for (int x=1; x<=MX; x++)
		for (int y=1; y<=MY; y++)
		for (int z=1; z<=MZ; z++) {
			theta_block+=phitot[(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)];
		}
		for (int x=1; x<=MX; x++)
		for (int y=1; y<=MY; y++)
		for (int z=1; z<=MZ; z++) {
			phitot[(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)]*=block[blocknr]/theta_block;
		        int length=MolMonList.size();
			for (int kk=0; kk<length; kk++)
			     phi[kk*M+(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)]*=block[blocknr]/theta_block;
		}
	}
}

void Molecule::SetThetaBlocks(int split) {
	int MX=lat->MX/split;
	int MY=lat->MY/split;
	int MZ=lat->MZ/split;
	int JX=lat->JX;
	int JY=lat->JY;
	Real theta_block;
	Real theta_tot=0;
	block.clear();
	for(int i=0; i<split; i++)
	for(int j=0; j<split; j++)
	for(int k=0; k<split; k++){
		theta_block=0;
		for (int x=1; x<=MX; x++)
		for (int y=1; y<=MY; y++)
		for (int z=1; z<=MZ; z++) {
			theta_block+=phitot[(i*MX+x)*JX+(j*MY+y)*JY+(k*MZ+z)];
		}
		theta_tot+=theta_block;
		block.push_back(theta_block);
	}
}

Real* Molecule::propagate_forward(Real* G1, int &s, int block, int generation, int M) {
NAMICS_DBG("1. propagate_forward for Mol " + name << endl);

	int N= n_mon[block];

	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];

		t=1;
		v0=t0=k0=0;

		if (s==first_s[generation]) {

			lat->Initiate(Gs+M,G1,Markov,M);
			lat->Initiate(Gs,G1,Markov,M); //not sure why this is done....
		} else {
			lat->propagate(Gs,G1,0,1,M); //assuming Gs contains previous end-point distribution on pos zero;
		}
		std::copy_n(Gs+M, M, Gg_f+n0*M); last_stored[block]=n0;
		s++;
		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {
			t++; s++;
			lat->propagate(Gs,G1,(k-1)%2,k%2,M);
			if (t>n) {
				t0++;
				if (t0 == n) t0 = ++v0;
				t = t0 + 1;
				k0 = k - t0 - 1;
			}

			if ((t == t0+1 && t0 == v0)
		  	 || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(k0+t0)))
		  	 || (2*(n-t+k) >= N-1)) {
				std::copy_n(Gs+(k%2)*M, M, Gg_f+(n0+t-1)*M);
				last_stored[block]=n0+t-1;
			}
		}
		if ((N)%2!=0) {
			std::copy_n(Gs+M, M, Gs);
		}
	} else {
		for (int k=0; k<N; k++) {
			if (s>first_s[generation]) {

				lat->propagate(Gg_f,G1,s-1,s,M);
			} else {
				lat->Initiate(Gg_f+first_s[generation]*M,G1,Markov,M);
			}
			 s++;
		}
	}
	if (save_memory) {
		return Gg_f+last_stored[block]*M;
	} else {
		 return Gg_f+(s-1)*M;
	}

}

void Molecule::propagate_backward(Real* G1, int &s, int block, int unity, int M) {
	(void)unity;
NAMICS_DBG("propagate_backward for Mol " + name << endl);

	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t,rk1;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];

		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {t++; if (t>n) { t0++; if (t0 == n) t0 = ++v0; t = t0 + 1; k0 = k - t0 - 1;}}
		for (k=N; k>=1; k--) {
			if (k==N) {
				if (s==chainlength-1) {
					std::copy_n(G1, M, Gg_b+(k%2)*M);
				} else {
					lat->propagate(Gg_b,G1,(k+1)%2,k%2,M);
				}
			} else {
				lat->propagate(Gg_b,G1,(k+1)%2,k%2,M);
			}
			t = k - k0;

			if (t == t0) {
				k0 += - n + t0;
				if (t0 == v0 ) {
					k0 -= ((n - t0)*(n - t0 + 1))/2;
				}
				t0 --;
				if (t0 < v0) {
					v0 = t0;
				}
				std::copy_n(Gg_f+(n0+t-1)*M, M, Gs+(t%2)*M);
				for (rk1=k0+t0+2; rk1<=k; rk1++) {
					t++;
					lat->propagate(Gs,G1,(t-1)%2,t%2,M);
					if (t == t0+1 || k0+n == k) {
						std::copy_n(Gs+(t%2)*M, M, Gg_f+(n0+t-1)*M);
					}
					if (t == n && k0+n < k) {
						t  = ++t0;
						k0 += n - t0;
					}
				}
				t = n;
			}
			lat->AddPhiS(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M,Gg_b+(k%2)*M,Markov,M);
			s--;
		}
		std::copy_n(Gg_b+M, M, Gg_b);
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) {
				lat->propagate(Gg_b,G1,(s+1)%2,s%2,M);
			} else {
				lat->Initiate(Gg_b+(s%2)*M,G1,Markov,M);
			}

			lat->AddPhiS(rho+molmon_nr[block]*M, Gg_f+(s*M), Gg_b+(s%2)*M,Markov, M);
			s--;
		}
	}
}



Real* Molecule::propagate_forward(Real* G1, int &s, int block, Real* P, int generation, int M) {
NAMICS_DBG("1. propagate_forward for Mol " + name << endl);

	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];
		if (s==first_s[generation]) {

			lat->Initiate(Gs+size*M,G1,Markov,M);
			//lat->Initiate(Gs,G1,Markov,M); //not necessary.
		} else {
			lat->propagateF(Gs,G1,P,0,1,M); //assuming Gs contains previous end-point distribution on pos zero;

		}
		s++;
		std::copy_n(Gs+M*size, M*size, Gg_f+n0*M*size);
		last_stored[block]=n0;

		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {
			t++; s++;
			lat->propagateF(Gs,G1,P,(k-1)%2,k%2,M);
			if (t>n) {
				t0++;
				if (t0 == n) { t0 = ++v0; cout <<"v0 >0 .... save memory may fail!" << endl; }
				t = t0 + 1;
				k0 = k - t0 - 1;
			}
			if ((t == t0+1 && t0 == v0)
		  	 || (t == t0+1 && ((n-t0)*(n-t0+1) >= N-1-2*(k0+t0)))
		  	 || (2*(n-t+k) >= N-1)) {
				std::copy_n(Gs+(k%2)*M*size, M*size, Gg_f+(n0+t-1)*M*size);
				last_stored[block]=n0+t-1;
			}
		}
		if ((N)%2!=0) {
			std::copy_n(Gs+M*size, M*size, Gs);
		}
	} else {
		for (int k=0; k<N; k++) {
			if (s>first_s[generation]) {
				lat->propagateF(Gg_f,G1,P,s-1,s,M);
			} else {
				lat->Initiate(Gg_f+first_s[generation]*M*size,G1,Markov,M);
			}
			 s++;
		}
	}
	if (save_memory) {
		return Gg_f+last_stored[block]*M*size;
	} else {
		 return Gg_f+(s-1)*M*size;
	}

}

void Molecule::propagate_backward(Real* G1, int &s, int block, Real* P, int& unity, int M) {
NAMICS_DBG("propagate_backward for Mol " + name << endl);
	int N= n_mon[block];
	if (save_memory) {
		int k,k0,t0,v0,t,rk1;
		int n=memory[block]; if (block>0) n-=memory[block-1];
		int n0=0; if (block>0) n0=memory[block-1];

		t=1;
		v0=t0=k0=0;
		for (k=2; k<=N; k++) {t++; if (t>n) { t0++; if (t0 == n) t0 = ++v0; t = t0 + 1; k0 = k - t0 - 1;}}
		for (k=N; k>=1; k--) {
			if (k==N) {
				if (s==chainlength-1) {
					lat->Initiate(Gg_b+(k%2)*M*size,G1,Markov,M);
				} else {
					if (unity==-1) {
						unity=0;
						Real* GB= (Real*) malloc(2*M*sizeof(Real));
						lat->Terminate(GB,Gg_b+((k+1)%2)*M*size,Markov,M);
						lat->propagate(GB,G1,0,1,M); //first step is freely joined
						lat->Initiate(Gg_b+(k%2)*M*size,GB+M,Markov,M);
						free(GB);
					} else {
						lat->propagateB(Gg_b,G1,P,(k+1)%2,k%2,M);
					}
				}
			} else {
				lat->propagateB(Gg_b,G1,P,(k+1)%2,k%2,M);
			}
			t = k - k0;
			if (t == t0) {
				k0 += - n + t0;
				if (t0 == v0 ) {
					k0 -= ((n - t0)*(n - t0 + 1))/2;
				}
				t0 --;
				if (t0 < v0) {
					v0 = t0;
				}
				std::copy_n(Gg_f+(n0+t-1)*M*size, M*size, Gs+(t%2)*M*size);
				for (rk1=k0+t0+2; rk1<=k; rk1++) {
					t++;
					lat->propagateF(Gs,G1,P,(t-1)%2,t%2,M);
					if (t == t0+1 || k0+n == k) {
						std::copy_n(Gs+(t%2)*M*size, M*size, Gg_f+(n0+t-1)*M*size);
					}
					if (t == n && k0+n < k) {
						t  = ++t0;
						k0 += n - t0;
					}
				}
				t = n;
			}

			lat->AddPhiS(rho+molmon_nr[block]*M,Gg_f+(n0+t-1)*M*size,Gg_b+(k%2)*M*size,Markov,M);
			s--;
		}
		std::copy_n(Gg_b+M*size, size*M, Gg_b);
	} else {
		for (int k=0; k<N; k++) {
			if (s<chainlength-1) {
				if (unity==-1) {
					unity=0;
					Real* GB= (Real*) malloc(2*M*sizeof(Real));
					lat->Terminate(GB,Gg_b+((s+1)%2)*M*size,Markov,M);
					lat->propagate(GB,G1,0,1,M); //first step is freely joined
					lat->Initiate(Gg_b+(s%2)*M*size,GB+M,Markov,M);
					free(GB);
				} else {
					lat->propagateB(Gg_b,G1,P,(s+1)%2,s%2,M);
				}
			} else {
				lat->Initiate(Gg_b+(s%2)*M*size,G1,Markov,M);
			}

			lat->AddPhiS(rho+molmon_nr[block]*M, Gg_f+s*M*size, Gg_b+(s%2)*M*size, Markov, M);
			s--;
		}
	}
}

bool Molecule::ComputePhi(Real* BETA,int id){
NAMICS_DBG("ComputePhi for Mol " + name << endl);
	bool success=true;
	int M=lat->M;
	if (id !=0) {
		int molmonlistlength= MolMonList.size();
		for (int i=0; i<molmonlistlength; i++)
		if (id==1) {
			for (int __i = 0; __i < (M); ++__i) (Seg[MolMonList[i]]->G1)[__i] = (Seg[MolMonList[i]]->G1)[__i] * (BETA)[__i];
		} else {
			for (int __i = 0; __i < (M); ++__i) (Seg[MolMonList[i]]->G1)[__i] = ((BETA)[__i] != 0) ? ((Seg[MolMonList[i]]->G1)[__i] / (BETA)[__i]) : 0;
		}
	}
	success=ComputePhi();


	if (id !=0) {
		int molmonlistlength=MolMonList.size();

		for (int i=0; i<molmonlistlength; i++) {
			if (id==1) {
				for (int __i = 0; __i < (M); ++__i) (phi+i*M)[__i] = ((BETA)[__i] != 0) ? ((phi+i*M)[__i] / (BETA)[__i]) : 0;
				for (int __i = 0; __i < (M); ++__i) (Seg[MolMonList[i]]->G1)[__i] = ((BETA)[__i] != 0) ? ((Seg[MolMonList[i]]->G1)[__i] / (BETA)[__i]) : 0;
			}  else {
				for (int __i = 0; __i < (M); ++__i) (phi+i*M)[__i] = (phi+i*M)[__i] * (BETA)[__i];
				for (int __i = 0; __i < (M); ++__i) (Seg[MolMonList[i]]->G1)[__i] = (Seg[MolMonList[i]]->G1)[__i] * (BETA)[__i];
			}
		}
	}
	return success;
}

bool Molecule::ComputePhi(){
NAMICS_DBG("ComputePhi for Molecule " + name << endl); //default computation for monomer only....
	int M=lat->M;
	bool success=true;
	std::copy_n(Seg[mon_nr[0]]->G1, M, phi);
	GN=lat->WeightedSum(phi);
	for (int __i = 0; __i < (M); ++__i) (phi)[__i] = (phi)[__i] * (Seg[mon_nr[0]]->G1)[__i];
	return success;
}


Real Molecule::fraction(int segnr){
NAMICS_DBG("fraction for mol_test " + name << endl); //default for monomer.
	int Nseg=0;
	int length = mon_nr.size();
	int i=0;
	if (ring) i++; //first segment is not counted in fraction;
	while (i<length) {
		if (segnr==mon_nr[i]) {Nseg+=n_mon[i];}
		i++;
	}
	return 1.0*Nseg/chainlength;
}
