#include "molecule.h"
#include "mol_branched.h"

#include <cctype>

namespace {

const std::vector<std::string>& MoleculeKeys() {
	static const std::vector<std::string> keys = {"freedom", "composition", "theta", "phibulk", "n", "B"};
	return keys;
}

bool ParseLinearCompositionLength(const std::string& composition, size_t& pos, int& chainlength) {
	chainlength = 0;
	while (pos < composition.size()) {
		if (composition[pos] != '(') return false;
		const size_t content_start = ++pos;
		int depth = 1;
		while (pos < composition.size() && depth > 0) {
			if (composition[pos] == '(') depth++;
			else if (composition[pos] == ')') depth--;
			pos++;
		}
		if (depth != 0 || pos >= composition.size() || !std::isdigit(static_cast<unsigned char>(composition[pos]))) return false;
		const size_t content_end = pos - 1;
		const size_t repeats_start = pos;
		while (pos < composition.size() && std::isdigit(static_cast<unsigned char>(composition[pos]))) pos++;
		const int repeats = ParseInt(composition.substr(repeats_start, pos - repeats_start), 0);
		if (repeats < 1) return false;

		const std::string content = composition.substr(content_start, content_end - content_start);
		int unit_length = 1;
		if (content.find('(') != std::string::npos) {
			size_t nested_pos = 0;
			if (!ParseLinearCompositionLength(content, nested_pos, unit_length) || nested_pos != content.size()) return false;
		}
		chainlength += unit_length * repeats;
	}
	return true;
}

bool IsMonomerComposition(const ParameterStore& parameters) {
	const std::string composition = parameters.value("composition", std::string{});
	if (composition.find('[') != std::string::npos || composition.find(']') != std::string::npos) return false;
	size_t pos = 0;
	int chainlength = 0;
	return ParseLinearCompositionLength(composition, pos, chainlength) && pos == composition.size() && chainlength == 1;
}

bool ExpandBrackets(const Input& input, const std::string& name, std::string& s) {
NAMICS_DBG("Molecule:: ExpandBrackets" << std::endl);
	bool success=true;
	if (s[0] != '(') {std::cout <<"illegal composition. Expects composition to start with a '(' in: " << s << std::endl; return false;}
	std::vector<int> open;
	std::vector<int> close;
	bool done=false; //now interpreted the (expanded) composition
	while (!done) { done = true;
		open.clear(); close.clear();
		if (!input.EvenBrackets(s,open,close)) {
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

bool Interpret(Molecule& mol, std::string s,int generation){
NAMICS_DBG("Molecule:: Interpret" << std::endl);
	if (s=="[") return true;
	bool success=true;
	std::vector<int>open;
	std::vector<int>close;
	mol.In->EvenBrackets(s,open,close);
	if (open.empty()) {
		std::cout << "In composition of mol '" + mol.name + "' an invalid token was found: " << s << std::endl;
		return false;
	}
	int k=0;
	int length=open.size();
	while (k<length) {
		std::string segname=s.substr(open[k]+1,close[k]-open[k]-1);
		int mnr=-1;
		int n_segments=mol.In->MonList.size();
		for (int i=0; i<n_segments; i++) {
			if (mol.Seg[i]->name ==segname) mnr=i;
		}
		if (mnr <0)  {std::cerr <<"In composition of mol '" + mol.name + "', segment name '" + segname + "' is not recognised"  << std::endl; success = false;
		throw "Composition Error";
		} else {

			int stored=mol.Gnr.size();
			if (stored>0) {//fragments at branchpoint need to be just 1 segment long.
				if (mol.Gnr[stored-1]<generation) {
					if (mol.n_mon[stored-1]>1) {
						mol.n_mon[stored-1]--;
						mol.n_mon.push_back(1);
						mol.mon_nr.push_back(mol.mon_nr[stored-1]);
						mol.Gnr.push_back(mol.Gnr[stored-1]);
					 	mol.last_b[mol.Gnr[stored-1]]++;
					}
				}
			}
			mol.mon_nr.push_back(mnr);
			mol.Gnr.push_back(generation);
			if (mol.first_s[generation] < 0) mol.first_s[generation]=mol.chainlength;
			if (mol.first_b[generation] < 0) mol.first_b[generation]=mol.mon_nr.size()-1;
			mol.last_b[generation]=mol.mon_nr.size()-1;
		}
		int nn = ParseInt(s.substr(close[k]+1,s.size()-close[k]-1),0);
		if (nn<1) {std::cout <<"In composition of mol '" + mol.name + "' the number of repeats should have values larger than unity " << std::endl; success=false; return success;
		} else {
			mol.n_mon.push_back(nn);
		}
		mol.chainlength +=nn; mol.last_s[generation]=mol.chainlength;
		k++;
	}
	return success;
}

bool GenerateTree(Molecule& mol, std::string s,int generation,int &pos, std::vector<int> open,std::vector<int> close) {
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
				mol.first_s.push_back(-1);
				mol.last_s.push_back(-1);
				mol.first_b.push_back(-1);
				mol.last_b.push_back(-1);
				success=GenerateTree(mol,s,new_generation,pos,open,close);
				if (!success) {std::cout <<"error in generate tree ." << std::endl; return success; }
				pos_close=pos_open+1;
			} else {
				pos=pos_close;
				success =Interpret(mol,ss,generation);
				if (!success)  {std::cout <<"error in interpret ." << std::endl; return success; }
			}
		} else {
			ss=s.substr(pos,pos_open-pos);
			pos=pos_open;
			success =Interpret(mol,ss,generation);
			if (!success) {std::cout <<"error in Interpret" << std::endl;  return success;}
			mol.first_s.push_back(-1);
			mol.last_s.push_back(-1);
			mol.first_b.push_back(-1);
			mol.last_b.push_back(-1);
			success=GenerateTree(mol,s,newgeneration,pos,open,close);
			if (!success) {std::cout <<"error in generate tree .." << std::endl; return success;}
		}
	}
	return success;
}

bool MakeMonList(Molecule& mol) {
NAMICS_DBG("Molecule:: MakeMonList" << std::endl);
	mol.MolMonList.clear();
	mol.molmon_nr.clear();
	bool success=true;
	int length = mol.mon_nr.size();
	int i=0;
	while (i<length) {
		if (!ContainsValue(mol.MolMonList,mol.mon_nr[i])) {
			if (mol.Seg[mol.mon_nr[i]]->freedom=="frozen") {
				success = false;
				std::cout << "In 'composition of mol " + mol.name + ", a segment was found with freedom 'frozen'. This is not permitted. " << std::endl;

			}
			mol.MolMonList.push_back(mol.mon_nr[i]);
		}
		i++;
	}
	i=0;
	int pos;
	while (i<length) {
		if (ContainsValue(mol.MolMonList,mol.mon_nr[i],&pos)) {mol.molmon_nr.push_back(pos);
		} else {std::cout <<"program error in mol PrepareForCalcualations" << std::endl; }
		i++;
	}
	return success;
}

bool Decomposition(Molecule& mol, std::string s){
NAMICS_DBG("Decomposition for Mol " + mol.name << std::endl);
	bool success = true;
	mol.chainlength=0;
	std::vector<int> open;
	std::vector<int> close;

	if (!ExpandBrackets(*mol.In, mol.name, s)) {success=false; return success;}
	if (!mol.In->EvenSquareBrackets(s,open,close)) {
		std::cout << "Error in composition of mol '" + mol.name + "'; the square brackets are not balanced in " << s << std::endl;
		success=false; return success;
	}
	const bool branched_molecule = !open.empty();
	int generation=0;
	int pos=0;
	mol.MolMonList.clear();
	int lopen;
	mol.Gnr.clear();
	mol.first_s.clear();
	mol.last_s.clear();
	mol.first_b.clear();
	mol.last_b.clear();
	mol.mon_nr.clear();
	mol.n_mon.clear();
	mol.molmon_nr.clear();
	mol.first_s.push_back(-1);
	mol.last_s.push_back(-1);
	mol.first_b.push_back(-1);
	mol.last_b.push_back(-1);
	if (branched_molecule) {
		lopen=open.size();
		for (int i=1; i<lopen; i++) {
			if ((open[i]-open[i-1])==1 || (close[i]-close[i-1])==1) {
				std::cout <<"In molecule " + mol.name + " in 'composition', two similar square brackets in a row '[[' or ']]' is not allowed" << std::endl;
				success=false;
				return success;
			}
		}
	}
	success = GenerateTree(mol,s,generation,pos,open,close);
	if (!success) {
		std::cout << "GenerateTree failed" << std::endl;
		return success;
	}
	if (branched_molecule) { //invert numbers;
		int g_length=mol.first_s.size();
		int length = mol.n_mon.size();
		int xxx;
		int ChainLength=mol.last_s[0];

		for (int i=0; i<g_length; i++) {
			mol.first_s[i] = ChainLength-mol.first_s[i]-1;
			mol.last_s[i] = ChainLength-mol.last_s[i]-1;
			xxx=mol.first_s[i]; mol.first_s[i]=mol.last_s[i]+1; mol.last_s[i]=xxx;
			mol.first_b[i]=length-mol.first_b[i]-1;
			mol.last_b[i]=length-mol.last_b[i]-1;
			xxx=mol.first_b[i]; mol.first_b[i]=mol.last_b[i]; mol.last_b[i]=xxx;
		}

		for (int i=0; i<length/2; i++) {
			xxx=mol.Gnr[i]; mol.Gnr[i]=mol.Gnr[length-1-i]; mol.Gnr[length-1-i]=xxx;
			xxx=mol.n_mon[i]; mol.n_mon[i]=mol.n_mon[length-1-i]; mol.n_mon[length-1-i]=xxx;
			xxx=mol.mon_nr[i]; mol.mon_nr[i]=mol.mon_nr[length-1-i]; mol.mon_nr[length-1-i]=xxx;
		}
		}

	success=MakeMonList(mol);
	return success;
}

bool CheckInput(Molecule& mol, int start_) {
NAMICS_DBG("Molecule:: CheckInput for mol " << mol.name << std::endl);
	mol.start=start_;
	mol.phibulk=0;
	mol.n=0;
	mol.theta=0;
	mol.norm=0;
NAMICS_DBG("CheckInput for Mol " + mol.name << std::endl);
	bool success=true;
	const auto& parameters = mol.In->Parameters("mol", mol.name, mol.start);
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(MoleculeKeys(), it.key())) continue;
		success = false;
		std::cout << "mol property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : MoleculeKeys()) std::cout << item << std::endl;
	}
	if (!success) return false;

	try {
		const std::string freedom_value = parameters.value("freedom", std::string{});
		const bool has_theta = parameters.contains("theta");
		const bool has_n = parameters.contains("n");
		const std::string composition_value = parameters.value("composition", std::string{});
		if (composition_value.size()==0) {
			std::cout << "For mol '" + mol.name + "' the definition of 'composition' is required" << std::endl;
			return false;
		}
		try {
			if (!Decomposition(mol, composition_value)) {
				std::cout << "For mol '" + mol.name + "' the composition is rejected. " << std::endl;
				return false;
			}
		} catch (const char* error) {
			std::cerr << error << std::endl;
			return false;
		}
		const bool pinned = mol.IsPinned();
		if (freedom_value.size()==0) {
			if (pinned) {
				std::cout <<"For mol " + mol.name + " the setting for 'freedom' was not set" << std::endl;
				return false;
			}
			std::cout <<"For mol " + mol.name + " the setting 'freedom' is expected: options: 'free' 'restricted' 'solvent' 'neutralizer' . Problem terminated " << std::endl;
			return false;
		}
		const bool allowed_freedom = freedom_value == "restricted" || (!pinned && (freedom_value == "free" || freedom_value == "solvent" || freedom_value == "neutralizer"));
		if (!allowed_freedom) {
			std::cout << "In mol " + mol.name + " the value for 'freedom' is not recognised " << std::endl;
			std::cout << "Select from: " << std::endl;
			if (!pinned) {
				std::cout << "free ; solvent ; neutralizer ; ";
			}
			std::cout << "restricted ; " << std::endl;
			return false;
		}
		mol.freedom = freedom_value;
		if (mol.freedom == "neutralizer" && !mol.IsCharged()) {
			std::cout << "Mol '" + mol.name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << std::endl;
			return false;
		}
		if (mol.freedom == "free") {
			if (!parameters.contains("phibulk")) {
				std::cout <<"In mol " + mol.name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<std::endl;
				return false;
			}
			mol.phibulk=parameters.at("phibulk").get<Real>();
			if (mol.phibulk < 0 || mol.phibulk >1) {
				std::cout << "In mol " + mol.name + ", the value of 'phibulk' is out of range 0 .. 1." << std::endl;
				return false;
			}
		}

		mol.B = 1;
		if (!pinned && parameters.contains("B")){
			mol.B=parameters.at("B").get<Real>();
			if (mol.B<1e-9) {
				std::cout <<"for Mol" + mol.name + " mobility B should have a posititve value. Default value B=1 is chosen. " << std::endl;
				mol.B=1;
			}
		}

		if (mol.freedom=="restricted") {
			if (!has_theta && !has_n) {
				std::cout <<"In mol " + mol.name + ", the setting 'freedom = restricted' should be combined with a value for 'theta' or 'n'; do not use both settings! "<<std::endl;
				return false;
			}
			if (has_theta && has_n) {
				std::cout <<"In mol " + mol.name + ", the setting 'freedom = restricted' does not allow both 'n' and 'theta' "<<std::endl;
				return false;
			}
			if (has_n) {mol.n=parameters.at("n").get<Real>(); mol.theta=mol.n*mol.chainlength;}
			if (has_theta) {mol.theta = parameters.at("theta").get<Real>(); mol.n=mol.theta/mol.chainlength;}
			if (mol.theta < 0 || (!pinned && mol.theta > mol.lat->volume)) {
				std::cout << "In mol " + mol.name + ", the value of 'n' or 'theta' is out of range." << std::endl;
				return false;
			}
		}
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in mol '" << mol.name << "': " << error.what() << std::endl;
		return false;
	}
	return true;
}

} // namespace

Molecule::Molecule(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_, std::string name_) {
	In=In_; Seg=Seg_; name=name_;
NAMICS_DBG("Constructor for Mol " + name << std::endl);
	lat=Lat_;
	all_molecule=false;
	B=1;

}

Molecule::~Molecule() {
	DeAllocateMemory();
}

void Molecule :: DeAllocateMemory(){
NAMICS_DBG("DeallocateMemory for Mol " + name << std::endl);
	if (!all_molecule) return;
	phi.clear();
	phi_ranked.clear();
	phitot.clear();
	Gg_f.clear();
	Gg_b.clear();
	UNITY.clear();
	all_molecule=false;
}

void Molecule:: AllocateMemory() {
	NAMICS_DBG("AllocateMemory in Mol " + name << std::endl);
	DeAllocateMemory();
	int M=lat->M;
	N=0;
	int length_ = mon_nr.size();
	for (int i=0; i<length_; i++) {N+=n_mon[i];}

	phi.assign(M * MolMonList.size(), 0);
	if (HasOutputProperty("phi_ranked")) phi_ranked.assign(M * N, 0);
	phitot.assign(M, 0);
	Gg_f.assign(M * N, 0);
	Gg_b.assign(2 * M, 0);
	UNITY.assign(M, 0);
	all_molecule=true;
}

bool Molecule:: PrepareForCalculations(std::span<const Real> KSAM) {
NAMICS_DBG("PrepareForCalculations in Mol " + name << std::endl);
	std::copy(KSAM.begin(), KSAM.end(), UNITY.begin());
	std::fill(phitot.begin(), phitot.end(), 0);
	std::fill(phi.begin(), phi.end(), 0);
	if (!phi_ranked.empty()) std::fill(phi_ranked.begin(), phi_ranked.end(), 0);
	return true;
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

namespace molecule_factory {

std::unique_ptr<Molecule> CreateChecked(const Input& input, Lattice* lat, std::span<const std::unique_ptr<Segment>> segments, const std::string& name, int start) {
	const auto& parameters = input.Parameters("mol", name, start);
	std::unique_ptr<Molecule> molecule;
	if (IsMonomerComposition(parameters)) molecule = std::make_unique<Molecule>(&input, lat, segments, name);
	else molecule = std::make_unique<mol_branched>(&input, lat, segments, name);
	if (!CheckInput(*molecule, start)) return nullptr;
	return molecule;
}

} // namespace molecule_factory

bool Molecule::HasOutputProperty(const std::string& property) const {
	const auto& problem = In->Start(start);
	if (!problem.is_object()) return false;
	const auto json_it = problem.find("json");
	if (json_it == problem.end() || !json_it->is_object()) return false;
	const auto mol_it = json_it->find("mol");
	if (mol_it == json_it->end() || !mol_it->is_object()) return false;
	const auto has_property = [&](ParameterStore::const_iterator it) {
		if (it == mol_it->end()) return false;
		const auto& value = it.value();
		if (value.is_string()) return value.get<std::string>() == property;
		if (value.is_array()) for (const auto& entry : value) if (entry.is_string() && entry.get<std::string>() == property) return true;
		return false;
	};
	return has_property(mol_it->find(name)) || has_property(mol_it->find("*"));
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
	if (!phi_ranked.empty()) OUTPUT["phi_ranked"] = {{"ranked_profile", 0}};
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
			lat->Initiate(Gg_f.data()+first_s[generation]*M,G1);
		}
		s++;
	}
	return Gg_f.data()+(s-1)*M;

}

void Molecule::propagate_backward(Real* G1, int &s, int block, int M) {
NAMICS_DBG("propagate_backward for Mol " + name << std::endl);

	int N= n_mon[block];
	for (int k=0; k<N; k++) {
		if (s<chainlength-1) {
			lat->propagate(Gg_b.data(),G1,(s+1)%2,s%2,M);
		} else {
			lat->Initiate(Gg_b.data()+(s%2)*M,G1);
		}

		lat->AddPhiS(phi.data()+molmon_nr[block]*M, Gg_f.data()+(s*M), Gg_b.data()+(s%2)*M);
		if (!phi_ranked.empty()) lat->AddPhiS(phi_ranked.data()+static_cast<size_t>(s)*M, Gg_f.data()+(s*M), Gg_b.data()+(s%2)*M);
		s--;
	}
}

bool Molecule::ComputePhi(){
NAMICS_DBG("ComputePhi for Molecule " + name << std::endl); // Default computation for a monomer.
	int M=lat->M;
	std::copy_n(Seg[mon_nr[0]]->G1.begin(), M, phi.begin());
	GN=lat->WeightedSum(phi.data());
	for (int __i = 0; __i < (M); ++__i) (phi)[__i] = (phi)[__i] * (Seg[mon_nr[0]]->G1)[__i];
	if (!phi_ranked.empty()) std::copy_n(phi.begin(), M, phi_ranked.begin());
	return true;
}


Real Molecule::fraction(int segnr){
NAMICS_DBG("fraction for Molecule " + name << std::endl); // Default for a monomer.
	int Nseg=0;
	int length = mon_nr.size();
	for (int i = 0; i < length; i++) {
		if (segnr==mon_nr[i]) {Nseg+=n_mon[i];}
	}
	return 1.0*Nseg/chainlength;
}
