#include "molecule.h"

Molecule::Molecule(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_, std::string name_) {
	In=In_; Seg=Seg_; name=name_;
NAMICS_DBG("Constructor for Mol " + name << std::endl);
	lat=Lat_;
	all_molecule=false;
	Markov =1;
	B=1;

}

Molecule::~Molecule() {
	DeAllocateMemory();
}

void Molecule :: DeAllocateMemory(){
NAMICS_DBG("DeallocateMemory for Mol " + name << std::endl);
	if (!all_molecule) return;
	phi.clear();
	phitot.clear();
	P.clear();
	Gg_f.clear();
	Gg_b.clear();
	UNITY.clear();
	all_molecule=false;
}

void Molecule:: AllocateMemory() {
NAMICS_DBG("AllocateMemory in Mol " + name << std::endl);
	DeAllocateMemory();
	int M=lat->M;
	if (Markov==2){// && lat->lattice_type == simple_cubic) {
		int FJC = lat->FJC;
		P.assign(FJC, 0);
		Real Q=0;
		KStiff=k_stiff;
		for (int k=0; k<FJC-1; k++) {
			P[k]=std::exp(-0.5*KStiff*(k*PIE/(FJC-1))*(k*PIE/(FJC-1)) ); //alternative to put U(theta)=-k(1-cos(theta))
			if (k>0) {
				if (lat->lattice_type==hexagonal) Q+= 2*P[k]; else Q+= 4*P[k]; //alternative is to use u_bend = Kstiff(1-cos(theta)), persistence length is l_p = b/ln <cos (theta)>
			} else Q=P[k];
		}
		P[FJC-1]=0; //Q+=P[FJC-1];
		if (lat->lattice_type==hexagonal&& !lat->stencil_full) Q*=2.0;
		for (int k=0; k<FJC; k++) { P[k]/=Q;
			std::cout << "P["<<k<<"] = " << P[k] << std::endl;
		}
	}

	N=0;
	int length_ = mon_nr.size();
	for (int i=0; i<length_; i++) {N+=n_mon[i];}

	phi.assign(M * MolMonList.size(), 0);
	phitot.assign(M, 0);
	Gg_f.assign(M * N * size, 0);
	Gg_b.assign(2 * M * size, 0);
	UNITY.assign(M * size, 0);
	all_molecule=true;
}

bool Molecule:: PrepareForCalculations(std::span<const Real> KSAM) {
NAMICS_DBG("PrepareForCalculations in Mol " + name << std::endl);
	std::copy(KSAM.begin(), KSAM.end(), UNITY.begin());
	bool success=true;
	std::fill(phitot.begin(), phitot.end(), 0);
	std::fill(phi.begin(), phi.end(), 0);


	return success;
}

bool Molecule::CheckInput(int start_, bool checking) {
NAMICS_DBG("Molecule:: CheckInput for mol " << name << std::endl);
	start=start_;
	phibulk=0;
	n=0;
	theta=0;
	norm=0;
NAMICS_DBG("CheckInput for Mol " + name << std::endl);
	bool success=true;
	const auto& parameters = In->Parameters("mol", name, start);
	static const std::vector<std::string> keys = {"freedom", "composition", "theta", "phibulk", "n", "Markov", "k_stiff", "B"};
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		success = false;
		std::cout << "mol property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : keys) std::cout << item << std::endl;
	}
	if (success) {
		try {
		const std::string freedom_value = parameters.value("freedom", std::string{});
		const bool has_theta = parameters.contains("theta");
		const bool has_n = parameters.contains("n");
		const std::string composition_value = parameters.value("composition", std::string{});
		if (composition_value.size()==0) {
			std::cout << "For mol '" + name + "' the definition of 'composition' is required" << std::endl;
			success = false;
		} else {
			try {
				if (!Decomposition(composition_value)) {
					std::cout << "For mol '" + name + "' the composition is rejected. " << std::endl;
					success=false;
				}
			} catch (const char* error) {
				std::cerr << error << std::endl;
				success = false;
			}
		}
		const bool pinned = IsPinned();
		if (checking) {
			if (freedom_value.size() > 0) freedom = freedom_value; else {
				std::cout <<"For molecule " << name << " no value for 'freedom' was found " << std::endl;
				success=false;
			}
			return success; //because moltype and freedom are known; as start <0 the checkinput can be terminated.
		}
		if (freedom_value.size()==0) {
			if (pinned) {
				std::cout <<"For mol " + name + " the setting for 'freedom' was not set" << std::endl;
				return false;
			}
			std::cout <<"For mol " + name + " the setting 'freedom' is expected: options: 'free' 'restricted' 'solvent' 'neutralizer' . Problem terminated " << std::endl;
			success = false;
		} else {
			const bool allowed_freedom = freedom_value == "restricted" || (!pinned && (freedom_value == "free" || freedom_value == "solvent" || freedom_value == "neutralizer"));
			if (!allowed_freedom) {
				std::cout << "In mol " + name + " the value for 'freedom' is not recognised " << std::endl;
				std::cout << "Select from: " << std::endl;
				if (!pinned) {
					std::cout << "free ; solvent ; neutralizer ; ";
				}
				std::cout << "restricted ; " << std::endl;
				if (pinned) return false;
				success=false;
			} else {
				freedom = freedom_value;
				if (freedom == "neutralizer" && !IsCharged()) {
					success=false;
					std::cout << "Mol '" + name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << std::endl;
				}
				if (freedom == "free") {
					if (!parameters.contains("phibulk")) {
						std::cout <<"In mol " + name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<std::endl;
						return false;
					}
					phibulk=parameters.at("phibulk").get<Real>();
					if (phibulk < 0 || phibulk >1) {
						std::cout << "In mol " + name + ", the value of 'phibulk' is out of range 0 .. 1." << std::endl;
						return false;
					}
				}

				if (!pinned) B=1;
				if (!pinned && parameters.contains("B")){
					B=parameters.at("B").get<Real>();
					if (B<1e-9) {
						std::cout <<"for Mol" + name + " mobility B should have a posititve value. Default value B=1 is chosen. " << std::endl;
						B=1;
					}
				}

				if (freedom=="restricted") {
					if (!has_theta && !has_n) {
						std::cout <<"In mol " + name + ", the setting 'freedom = restricted' should be combined with a value for 'theta' or 'n'; do not use both settings! "<<std::endl;
						success=false;
					} else if (has_theta && has_n) {
						std::cout <<"In mol " + name + ", the setting 'freedom = restricted' does not allow both 'n' and 'theta' "<<std::endl;
						success=false;
					} else {
						if (has_n) {n=parameters.at("n").get<Real>();theta=n*chainlength;}
						if (has_theta) {theta = parameters.at("theta").get<Real>();n=theta/chainlength;}
						if (theta < 0 || (!pinned && theta > lat->volume)) {
							std::cout << "In mol " + name + ", the value of 'n' or 'theta' is out of range." << std::endl;
							success=false;
						}
					}
				}
			}
		}
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid json type in mol '" << name << "': " << error.what() << std::endl;
			success = false;
		}
	}
	Markov=1;
	if (parameters.contains("Markov")) Markov=parameters.at("Markov").get<int>();
	if (Markov<1 || Markov>2) {
		std::cout <<" Integer value for 'Markov' is by default 1 and may be set to 2 for some mol_types and fjc-choices only. Markov value out of bounds. Proceed with caution. " << std::endl;
		success = false;
	}
	if (Markov==2) lat->Markov=2;
	k_stiff=lat->k_stiff; //pick up 'default' value from lattice.
	if (parameters.contains("k_stiff")) {
		k_stiff=parameters.at("k_stiff").get<Real>();
		if (k_stiff<0 || k_stiff>10) {
			success =false;
			std::cout <<" Real value for 'k_stiff' out of bounds (0 < k_stiff < 10). " << std::endl;
			std::cout <<" For Markov == 2: u_bend (theta) = 0.5 k_stiff theta^2, where 'theta' is angle for bond direction deviating from the straight direction. " <<std::endl;
			std::cout <<" You may interpret 'k_stiff' as the molecular 'persistence length' " << std::endl;
			std::cout <<" k_stiff is a 'default value'. Use molecular specific values to overrule the default when appropriate (future implementation....) " << std::endl;
		}
		if (lat->fjc>1 && lat->gradients>1) {
			success=false;
			std::cout <<" Work in progress.... Currently, Markov == 2 is implemented in gradients>1 for FJC_choices = 3 " << std::endl;
		}
	}
	size=0;
	if (Markov ==2) {
		if (lat->gradients==1) size = lat->FJC;
		if (lat->gradients==2) { //assume fjc=1...
		    if (lat->stencil_full) {
				lat->stencil_full=false; std::cout <<"Warning: stencil_full is set to false" << std::endl;
			}
			if (lat->lattice_type==hexagonal ) size =  12;
			if (lat->lattice_type==simple_cubic ) size = 2*lat->FJC-1;
			if (lat->lattice_type==hexagonal) { success=false;
				std::cout <<"Warning: Markov = 2 & gradients=2 lattice_type = hexagonal...not certified. Caution recommended, even without obvious error messages..." << std::endl;
			}
		}
		if (lat->gradients==3) {
		    if (lat->stencil_full) {
				lat->stencil_full=false; std::cout <<"Warning: stencil_full is set to false" << std::endl;
			}

			if (lat->lattice_type==hexagonal) size=12;
			if (lat->lattice_type==simple_cubic) size = 6;
			if (lat->lattice_type==hexagonal) { success=false;
				std::cout <<"Warning: Markov = 2 & gradients=3 lattice_type = hexagonal...not certified!!!. Caution recommended, even without obvious error messages...." << std::endl;
			}
		}
	} else size = 1;
	if (size==0) {success=false; std::cout <<"Attention: size in molecule is not set; combination gradients (1,2,3),  Markov=2, stencil_full (true,false), lattice_type (hexagonal, simple_cubic) not implemented" << std::endl;}

	return success;
}

bool Molecule::ExpandBrackets(std::string &s) {
NAMICS_DBG("Molecule:: ExpandBrackets" << std::endl);
	bool success=true;
	if (s[0] != '(') {std::cout <<"illegal composition. Expects composition to start with a '(' in: " << s << std::endl; return false;}
	std::vector<int> open;
	std::vector<int> close;
	bool done=false; //now interpreted the (expanded) composition
	while (!done) { done = true;
		open.clear(); close.clear();
		if (!In->EvenBrackets(s,open,close)) {
			std::cout << "s : " << s << std::endl;
			std::cout << "In composition of mol '" + name + "' the backets are not balanced."<<std::endl; success=false; return success;
		 }
		int length=open.size();
		int pos_open;
		int pos_close;
		int pos_low=0;
		int i_open=0; {pos_open=open[0]; pos_low=open[0];}
		int i_close=0; pos_close=close[0];
		if (pos_open > pos_close) {std::cout << "Brackets open in composition not correct" << std::endl; return false;}
		while (i_open < length-1 && done) {
			i_open++;
			pos_open=open[i_open];
			if (pos_open < pos_close && done) {
				i_close++; if (i_close<length) pos_close=close[i_close];
			} else {
				if (pos_low==open[i_open-1] ) {
					pos_low=pos_open;
					i_close++; if (i_close<length) pos_close=close[i_close];
					if (pos_open > pos_close) {std::cout << "Brackets open in composition not correct" << std::endl; return false;}
				} else {
					done=false;

					int x=ParseInt(s.substr(pos_close+1),-1);
						if (x<1) {
								std::cout <<"Number of repeats must be a positive integer in composition at pos : " << pos_close+1 << " for: " << s << std::endl; return false;
						}
					std::string sA,sB,sC;
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
			std::string sA,sB,sC;
			sA=s.substr(0,pos_low);
			sB=s.substr(pos_low+1,pos_close-pos_low-1);
			sC="";
			s=sA;for (int k=0; k<x; k++) s.append(sB); s.append(sC);
		}

	}
	if (s[s.size()-1]==']') {std::cout <<"illegal composition. Composition can not end with a ']' in: " << s << std::endl; return false;}
	return success;
}

bool Molecule::Interpret(std::string s,int generation){
NAMICS_DBG("Molecule:: Interpret" << std::endl);
	if (s=="[") return true;
	bool success=true;
	std::vector<int>open;
	std::vector<int>close;
	In->EvenBrackets(s,open,close);
	if (open.empty()) {
		std::cout << "In composition of mol '" + name + "' an invalid token was found: " << s << std::endl;
		return false;
	}
	int k=0;
	int length=open.size();
	while (k<length) {
		std::string segname=s.substr(open[k]+1,close[k]-open[k]-1);
		int mnr=-1;
		int n_segments=In->MonList.size();
		for (int i=0; i<n_segments; i++) {
			if (Seg[i]->name ==segname) mnr=i;
		}
		if (mnr <0)  {std::cerr <<"In composition of mol '" + name + "', segment name '" + segname + "' is not recognised"  << std::endl; success = false;
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
		if (nn<1) {std::cout <<"In composition of mol '" + name + "' the number of repeats should have values larger than unity " << std::endl; success=false; return success;
		} else {
			n_mon.push_back(nn);
		}
		chainlength +=nn; last_s[generation]=chainlength;
		k++;
	}
	return success;
}

bool Molecule::GenerateTree(std::string s,int generation,int &pos, std::vector<int> open,std::vector<int> close) {
NAMICS_DBG("Molecule:: GenerateTree" << std::endl);
	bool success=true;
	std::string ss;
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
				if (!success) {std::cout <<"error in generate tree ." << std::endl; return success; }
				pos_close=pos_open+1;
			} else {
				pos=pos_close;
				success =Interpret(ss,generation);
				if (!success)  {std::cout <<"error in interpret ." << std::endl; return success; }
			}
		} else {
			ss=s.substr(pos,pos_open-pos);
			pos=pos_open;
			success =Interpret(ss,generation);
			if (!success) {std::cout <<"error in Interpret" << std::endl;  return success;}
			first_s.push_back(-1);
			last_s.push_back(-1);
			first_b.push_back(-1);
			last_b.push_back(-1);
			success=GenerateTree(s,newgeneration,pos,open,close);
			if (!success) {std::cout <<"error in generate tree .." << std::endl; return success;}
		}
	}
	return success;
}

bool Molecule::Decomposition(std::string s){
NAMICS_DBG("Decomposition for Mol " + name << std::endl);
	bool success = true;
	MolType=linear;
	chainlength=0;
	std::vector<int> open;
	std::vector<int> close;

	if (!ExpandBrackets(s)) {success=false; return success;}
	if (!In->EvenSquareBrackets(s,open,close)) {
		std::cout << "Error in composition of mol '" + name + "'; the square brackets are not balanced in " << s << std::endl;
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
				std::cout <<"In molecule " + name + " in 'composition', two similar square brackets in a row '[[' or ']]' is not allowed" << std::endl;
				success=false;
				return success;
			}
		}
	}
	success = GenerateTree(s,generation,pos,open,close);
	if (!success) {
		std::cout << "GenerateTree failed" << std::endl;
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
		}

	success=MakeMonList();
	if (chainlength==1) MolType=monomer;
	return success;
}

bool Molecule:: MakeMonList(void) {
NAMICS_DBG("Molecule:: MakeMonList" << std::endl);
	MolMonList.clear();
	bool success=true;
	int length = mon_nr.size();
	int i=0;
	while (i<length) {
		if (!ContainsValue(MolMonList,mon_nr[i])) {
			if (Seg[mon_nr[i]]->freedom=="frozen") {
				success = false;
				std::cout << "In 'composition of mol " + name + ", a segment was found with freedom 'frozen'. This is not permitted. " << std::endl;

			}
			MolMonList.push_back(mon_nr[i]);
		}
		i++;
	}
	i=0;
	int pos;
	while (i<length) {
		if (ContainsValue(MolMonList,mon_nr[i],&pos)) {molmon_nr.push_back(pos);
		} else {std::cout <<"program error in mol PrepareForCalcualations" << std::endl; }
		i++;
	}
	return success;
}

bool Molecule::IsPinned() {
NAMICS_DBG("IsPinned for Mol " + name << std::endl);
	int length=MolMonList.size();
	for (int i=0; i<length; i++) {
		if (Seg[MolMonList[i]]->freedom=="pinned") return true;
	}
	return false;
}

Real Molecule::Charge() {
NAMICS_DBG("Molecule:: Charge" << std::endl);
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
NAMICS_DBG("IsCharged for Mol " + name << std::endl);
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

void Molecule::PushOutput() {
NAMICS_DBG("PushOutput for Mol " + name << std::endl);
	const auto& parameters = In->Parameters("mol", name, start);
	OUTPUT = nlohmann::ordered_json::object();
	OUTPUT["composition"] = parameters.value("composition", std::string{});
	OUTPUT["freedom"] = freedom;
	if (freedom == "free") theta = lat->WeightedSum(phitot.data());
	OUTPUT["Markov"] = Markov;
	OUTPUT["k_stiff"] = k_stiff;
	if (lat->gradients==3) {
		int MZ=lat->MZ;
		int MY=lat->MY;
		int MX=lat->MX;
		int JX=lat->JX;
		int JY=lat->JY;
		for (int z = 1; z < MZ + 1 && z <= 20; z++) {
			Real phiz=0;
			for (int x=1; x<MX+1; x++) for (int y=1;y<MY+1;y++) {
				phiz +=phitot[x*JX+y*JY+z];
			}
			phiz /= MX*MY;
			OUTPUT["phiz[" + std::to_string(z) + "]"] = phiz;
		}
	}
	if (Markov == 2) {
		for (int k = 0; k < size && k < 5; k++) {
			OUTPUT["P[" + std::to_string(k) + "]"] = P[k];
		}
	}
	OUTPUT["Rg"] = std::pow((lat->Moment(phitot.data(),0.0,2) / chainlength), 0.5);
	lat->remove_bounds(phitot.data());
	theta=lat->WeightedSum(phitot.data());
	OUTPUT["theta"] = theta;
	Real thetaexc=theta-lat->volume*phibulk;
	OUTPUT["theta_exc"] = thetaexc;
	OUTPUT["n_exc"] = thetaexc / chainlength;
	OUTPUT["nexc"] = thetaexc / chainlength;
	OUTPUT["thetaexc"] = thetaexc;
	OUTPUT["n"] = n;
	OUTPUT["chainlength"] = chainlength;
	OUTPUT["phibulk"] = phibulk;
	OUTPUT["Mu"] = Mu;
	OUTPUT["mu"] = Mu;
	OUTPUT["MU"] = Mu;
	if (lat->gradients==3) {
		Real TrueVolume=lat->MX*lat->MY*lat->MZ;
		Real Volume_particles=0;
		int num_of_seg=In->MonList.size();
		for (int i=0; i<num_of_seg; i++) {
			if (Seg[i]->freedom=="frozen") {
				for (int __j = 0; __j < lat->M; ++__j) Volume_particles += Seg[i]->MASK[__j];
			}
		}
		OUTPUT["Gamma"] = theta - (TrueVolume - Volume_particles) * phibulk;
	}
	if (chainlength==1) {
		int seg=MolMonList[0];
		if (Seg[seg]->ns >1) {
			if (mu_state.size() ==0) for (int i=0; i<Seg[seg]->ns; i++) mu_state.push_back(Mu);
			for (int i=0; i<Seg[seg]->ns; i++) {
				mu_state[i]+=std::log(Seg[seg]->state_alphabulk[i]);
				OUTPUT["mu-" + Seg[seg]->state_name[i]] = mu_state[i];
			}
		}
	}
	int M=lat->M;
	Real phimax=phitot[M/2];
	bool maxfound=false;
	int i=M/2;
	while (!maxfound) {
		i++;
		if (phitot[i]> phimax ) phimax =phitot[i]; else maxfound=true;
	}
	maxfound=false;
	i=M/2;
	while (!maxfound) {
		i--;
		if (phitot[i]> phimax ) phimax =phitot[i]; else maxfound=true;
	}

	OUTPUT["phiMax"] = phimax;
	OUTPUT["GN"] = GN;
	OUTPUT["norm"] = norm;
	OUTPUT["phi"] = {{"profile", 0}};
	for (size_t i = 0; i < MolMonList.size(); i++) {
		OUTPUT["phi_" + Seg[MolMonList[i]]->name] = {{"profile", static_cast<int>(i) + 1}};
	}
}

std::span<Real> Molecule::GetPointer(int profile) {
NAMICS_DBG("GetPointer for Mol " + name << std::endl);
	const int M = lat->M;
	if (profile == 0) {
		lat->set_bounds(phitot.data());
		return phitot;
	}
	if (profile > 0 && profile <= static_cast<int>(MolMonList.size())) {
		auto data = std::span<Real>(phi).subspan(static_cast<size_t>(profile - 1) * M, static_cast<size_t>(M));
		lat->set_bounds(data.data());
		return data;
	}
	return {};
}

Real* Molecule::propagate_forward(Real* G1, int &s, int block, int generation, int M) {
NAMICS_DBG("1. propagate_forward for Mol " + name << std::endl);

	int N= n_mon[block];
	for (int k=0; k<N; k++) {
		if (s>first_s[generation]) {
			lat->propagate(Gg_f.data(),G1,s-1,s,M);
		} else {
			lat->Initiate(Gg_f.data()+first_s[generation]*M,G1,Markov,M);
		}
		s++;
	}
	return Gg_f.data()+(s-1)*M;

}

void Molecule::propagate_backward(Real* G1, int &s, int block, int unity, int M) {
	(void)unity;
NAMICS_DBG("propagate_backward for Mol " + name << std::endl);

	int N= n_mon[block];
	for (int k=0; k<N; k++) {
		if (s<chainlength-1) {
			lat->propagate(Gg_b.data(),G1,(s+1)%2,s%2,M);
		} else {
			lat->Initiate(Gg_b.data()+(s%2)*M,G1,Markov,M);
		}

		lat->AddPhiS(phi.data()+molmon_nr[block]*M, Gg_f.data()+(s*M), Gg_b.data()+(s%2)*M,Markov, M);
		s--;
	}
}



Real* Molecule::propagate_forward(Real* G1, int &s, int block, Real* P, int generation, int M) {
NAMICS_DBG("1. propagate_forward for Mol " + name << std::endl);

	int N= n_mon[block];
	for (int k=0; k<N; k++) {
		if (s>first_s[generation]) {
			lat->propagateF(Gg_f.data(),G1,P,s-1,s,M);
		} else {
			lat->Initiate(Gg_f.data()+first_s[generation]*M*size,G1,Markov,M);
		}
		s++;
	}
	return Gg_f.data()+(s-1)*M*size;

}

void Molecule::propagate_backward(Real* G1, int &s, int block, Real* P, int& unity, int M) {
NAMICS_DBG("propagate_backward for Mol " + name << std::endl);
	int N= n_mon[block];
	for (int k=0; k<N; k++) {
		if (s<chainlength-1) {
			if (unity==-1) {
				unity=0;
				std::vector<Real> GB(2 * M, 0);
				lat->Terminate(GB.data(),Gg_b.data()+((s+1)%2)*M*size,Markov,M);
				lat->propagate(GB.data(),G1,0,1,M); //first step is freely joined
				lat->Initiate(Gg_b.data()+(s%2)*M*size,GB.data()+M,Markov,M);
			} else {
				lat->propagateB(Gg_b.data(),G1,P,(s+1)%2,s%2,M);
			}
		} else {
			lat->Initiate(Gg_b.data()+(s%2)*M*size,G1,Markov,M);
		}

		lat->AddPhiS(phi.data()+molmon_nr[block]*M, Gg_f.data()+s*M*size, Gg_b.data()+(s%2)*M*size, Markov, M);
		s--;
	}
}

bool Molecule::ComputePhi(){
NAMICS_DBG("ComputePhi for Molecule " + name << std::endl); //default computation for monomer only....
	int M=lat->M;
	bool success=true;
	std::copy_n(Seg[mon_nr[0]]->G1.begin(), M, phi.begin());
	GN=lat->WeightedSum(phi.data());
	for (int __i = 0; __i < (M); ++__i) (phi)[__i] = (phi)[__i] * (Seg[mon_nr[0]]->G1)[__i];
	return success;
}


Real Molecule::fraction(int segnr){
NAMICS_DBG("fraction for mol_test " + name << std::endl); //default for monomer.
	int Nseg=0;
	int length = mon_nr.size();
	for (int i = 0; i < length; i++) {
		if (segnr==mon_nr[i]) {Nseg+=n_mon[i];}
	}
	return 1.0*Nseg/chainlength;
}
