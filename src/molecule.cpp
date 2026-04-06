#include "molecule.h"

#include <cctype>

namespace {

const std::vector<std::string>& MoleculeKeys() {
	static const std::vector<std::string> keys = {"freedom", "composition", "theta", "phibulk", "n", "B"};
	return keys;
}

using Topology = Molecule::Topology;
using Node = Molecule::Node;

void SkipWhitespace(const std::string& s, size_t& pos) {
	while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
}

bool ParseRepeat(const std::string& s, size_t& pos, int& repeats) {
	SkipWhitespace(s, pos);
	if (pos >= s.size() || !std::isdigit(static_cast<unsigned char>(s[pos]))) return false;
	const size_t start = pos;
	while (pos < s.size() && std::isdigit(static_cast<unsigned char>(s[pos]))) ++pos;
	repeats = ParseInt(s.substr(start, pos - start), 0);
	return repeats > 0;
}

int SegmentIndex(const std::unordered_map<std::string, int>& lookup, const std::string& name) {
	const auto it = lookup.find(name);
	return it != lookup.end() ? it->second : -1;
}

bool ParseChain(const Molecule& mol, const std::string& s, size_t& pos, char terminator, Topology& chain, const std::unordered_map<std::string, int>& lookup, bool& branched);

bool ParseRawSegment(const Molecule& mol, const std::string& s, size_t& pos, Node& node, const std::unordered_map<std::string, int>& lookup) {
	const size_t start = pos;
	while (pos < s.size() && s[pos] != ')') {
		if (s[pos] == '(' || s[pos] == '[' || s[pos] == ']') return false;
		++pos;
	}
	if (pos >= s.size() || s[pos] != ')') return false;
	const std::string segname = s.substr(start, pos - start);
	++pos;
	const int seg_index = SegmentIndex(lookup, segname);
	if (seg_index < 0) {
		std::cerr << "In composition of mol '" + mol.name + "', segment name '" + segname + "' is not recognised" << std::endl;
		return false;
	}
	node.segment = seg_index;
	return true;
}

bool ParseChain(const Molecule& mol, const std::string& s, size_t& pos, char terminator, Topology& chain, const std::unordered_map<std::string, int>& lookup, bool& branched) {
	while (pos < s.size()) {
		SkipWhitespace(s, pos);
		if (pos >= s.size()) break;
		if (terminator != '\0' && s[pos] == terminator) {
			++pos;
			return true;
		}
		if (s[pos] == '[') {
			if (chain.empty()) {
				std::cout << "In composition of mol '" + mol.name + "' a branch was found before any segment." << std::endl;
				return false;
			}
			branched = true;
			++pos;
			Topology branch;
			if (!ParseChain(mol, s, pos, ']', branch, lookup, branched)) {
				std::cout << "In composition of mol '" + mol.name + "' the branch brackets are not balanced." << std::endl;
				return false;
			}
			if (branch.empty()) {
				std::cout << "In composition of mol '" + mol.name + "' contains an empty branch." << std::endl;
				return false;
			}
			chain.back().branches.push_back(std::move(branch));
			continue;
		}
		if (s[pos] != '(') {
			if (terminator == '\0') {
				std::cout << "In composition of mol '" + mol.name + "' an invalid token was found: " << s.substr(pos) << std::endl;
			}
			return false;
		}
		++pos;
		SkipWhitespace(s, pos);
		if (pos >= s.size()) return false;
		Node node;
		if (s[pos] == '(' || s[pos] == '[') {
			Topology group;
			if (!ParseChain(mol, s, pos, ')', group, lookup, branched)) {
				std::cout << "In composition of mol '" + mol.name + "' the brackets are not balanced." << std::endl;
				return false;
			}
			if (group.empty()) {
				std::cout << "In composition of mol '" + mol.name + "' contains an empty group." << std::endl;
				return false;
			}
			int repeats = 0;
			if (!ParseRepeat(s, pos, repeats)) {
				std::cout << "In composition of mol '" + mol.name + "' the number of repeats should have values larger than unity " << std::endl;
				return false;
			}
			chain.reserve(chain.size() + group.size() * static_cast<size_t>(repeats));
			for (int i = 0; i < repeats; ++i) chain.insert(chain.end(), group.begin(), group.end());
			continue;
		}
		if (!ParseRawSegment(mol, s, pos, node, lookup)) {
			std::cout << "In composition of mol '" + mol.name + "' an invalid token was found: " << s << std::endl;
			return false;
		}
		int repeats = 0;
		if (!ParseRepeat(s, pos, repeats)) {
			std::cout << "In composition of mol '" + mol.name + "' the number of repeats should have values larger than unity " << std::endl;
			return false;
		}
		node.repeat = repeats;
		chain.push_back(std::move(node));
	}
	return terminator == '\0';
}

bool ParseTopology(Molecule& mol, const std::string& composition, bool& branched) {
	mol.topology.clear();
	mol.MolMonList.clear();
	mol.molmon_nr.clear();
	mol.Gnr.clear();
	mol.first_s.clear();
	mol.last_s.clear();
	mol.first_b.clear();
	mol.last_b.clear();
	mol.mon_nr.clear();
	mol.n_mon.clear();
	mol.chainlength = 0;
	branched = false;

	std::unordered_map<std::string, int> lookup;
	for (size_t i = 0; i < mol.Seg.size(); ++i) lookup[mol.Seg[i]->name] = static_cast<int>(i);

	size_t pos = 0;
	if (!ParseChain(mol, composition, pos, '\0', mol.topology, lookup, branched)) return false;
	SkipWhitespace(composition, pos);
	if (pos != composition.size() || mol.topology.empty()) {
		std::cout << "In composition of mol '" + mol.name + "' an invalid token was found: " << composition.substr(pos) << std::endl;
		return false;
	}
	return true;
}

void EnsureGeneration(Molecule& mol, int generation) {
	if (generation < static_cast<int>(mol.first_s.size())) return;
	mol.first_s.resize(generation + 1, -1);
	mol.last_s.resize(generation + 1, -1);
	mol.first_b.resize(generation + 1, -1);
	mol.last_b.resize(generation + 1, -1);
}

void BuildTopologyCaches(Molecule& mol, const Topology& chain, int generation, int& next_generation) {
	EnsureGeneration(mol, generation);
	for (const auto& node : chain) {
		const int block_index = static_cast<int>(mol.mon_nr.size());
		mol.mon_nr.push_back(node.segment);
		mol.n_mon.push_back(node.repeat);
		mol.Gnr.push_back(generation);
		if (mol.first_s[generation] < 0) mol.first_s[generation] = mol.chainlength;
		if (mol.first_b[generation] < 0) mol.first_b[generation] = block_index;
		mol.last_b[generation] = block_index;
		mol.chainlength += node.repeat;
		mol.last_s[generation] = mol.chainlength;
		for (const auto& branch : node.branches) {
			const int branch_generation = next_generation++;
			BuildTopologyCaches(mol, branch, branch_generation, next_generation);
		}
	}
}

bool MakeMonList(Molecule& mol) {
NAMICS_DBG("Molecule:: MakeMonList" << std::endl);
	mol.MolMonList.clear();
	mol.molmon_nr.clear();
	mol.MolMonList.reserve(mol.Seg.size());
	bool success = true;
	std::vector<int> segment_positions(mol.Seg.size(), -1);
	for (int seg : mol.mon_nr) {
		if (segment_positions[seg] >= 0) continue;
		if (mol.Seg[seg]->freedom == "frozen") {
			success = false;
			std::cout << "In 'composition of mol " + mol.name + ", a segment was found with freedom 'frozen'. This is not permitted. " << std::endl;
		}
		segment_positions[seg] = static_cast<int>(mol.MolMonList.size());
		mol.MolMonList.push_back(seg);
	}
	mol.molmon_nr.reserve(mol.mon_nr.size());
	for (int seg : mol.mon_nr) {
		const int pos = segment_positions[seg];
		if (pos >= 0) {
			mol.molmon_nr.push_back(pos);
		} else {
			std::cout <<"program error in mol PrepareForCalcualations" << std::endl;
		}
	}
	return success;
}

bool Decomposition(Molecule& mol, std::string s){
NAMICS_DBG("Decomposition for Mol " + mol.name << std::endl);
	bool branched = false;
	if (!ParseTopology(mol, s, branched)) return false;
	int next_generation = 1;
	BuildTopologyCaches(mol, mol.topology, 0, next_generation);
	if (branched) {
		const int g_length = static_cast<int>(mol.first_s.size());
		const int length = static_cast<int>(mol.n_mon.size());
		int xxx;
		const int ChainLength = mol.last_s[0];

		for (int i = 0; i < g_length; i++) {
			mol.first_s[i] = ChainLength - mol.first_s[i] - 1;
			mol.last_s[i] = ChainLength - mol.last_s[i] - 1;
			xxx = mol.first_s[i]; mol.first_s[i] = mol.last_s[i] + 1; mol.last_s[i] = xxx;
			mol.first_b[i] = length - mol.first_b[i] - 1;
			mol.last_b[i] = length - mol.last_b[i] - 1;
			xxx = mol.first_b[i]; mol.first_b[i] = mol.last_b[i]; mol.last_b[i] = xxx;
		}

		for (int i = 0; i < length / 2; i++) {
			xxx = mol.Gnr[i]; mol.Gnr[i] = mol.Gnr[length - 1 - i]; mol.Gnr[length - 1 - i] = xxx;
			xxx = mol.n_mon[i]; mol.n_mon[i] = mol.n_mon[length - 1 - i]; mol.n_mon[length - 1 - i] = xxx;
			xxx = mol.mon_nr[i]; mol.mon_nr[i] = mol.mon_nr[length - 1 - i]; mol.mon_nr[length - 1 - i] = xxx;
		}
	}

	return MakeMonList(mol);
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
	for (int count : n_mon) N += count;

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
	return true;
}

void Molecule::FinalizeOutputs() {
	if (phi_ranked.empty()) return;
	std::fill(phi_ranked.begin(), phi_ranked.end(), 0);
	ComputePhiRanked(std::span<Real>(phi_ranked));
	int M = lat->M;
	const size_t Msize = M;
	size_t offset = 0;
	for (size_t b = 0; b < mon_nr.size(); ++b) {
		const auto g1 = std::span<const Real>(Seg[mon_nr[b]]->G1);
		for (int repeat = 0; repeat < n_mon[b]; ++repeat, offset += Msize) {
			auto phi = std::span<Real>(phi_ranked).subspan(offset, Msize);
			for (int i = 0; i < M; ++i) phi[i] = g1[i] != 0 ? phi[i] / g1[i] : 0;
			if (norm > 0) {
				for (int i = 0; i < M; ++i) phi[i] *= norm;
			}
		}
	}

	const auto monomers = SegmentIndices();
	for (size_t k = 0; k < monomers.size(); ++k) {
		const int seg = monomers[k];
		auto phi_block = std::span<Real>(phi).subspan(k * Msize, Msize);
		const auto g1 = std::span<const Real>(Seg[seg]->G1);
		for (int i = 0; i < M; ++i) phi_block[i] = g1[i] != 0 ? phi_block[i] / g1[i] : 0;
		if (norm > 0) {
			for (int i = 0; i < M; ++i) phi_block[i] *= norm;
		}
	}
}

bool Molecule::IsPinned() {
NAMICS_DBG("IsPinned for Mol " + name << std::endl);
	for (int seg : SegmentIndices()) {
		if (Seg[seg]->freedom == "pinned") return true;
	}
	return false;
}

Real Molecule::Charge() {
NAMICS_DBG("Molecule:: Charge" << std::endl);
	Real charge=0;
	ForEachNode([&](const Node& node, int) {
		const auto& seg = *Seg[node.segment];
		if (seg.state_name.size() > 1) {
			for (size_t i = 0; i < seg.state_name.size(); ++i) charge += seg.state_alphabulk[i] * seg.state_valence[i] * node.repeat;
		} else {
			charge += seg.valence * node.repeat;
		}
	});
	return charge/chainlength;
}

namespace molecule_factory {

std::unique_ptr<Molecule> CreateChecked(const Input& input, Lattice* lat, std::span<const std::unique_ptr<Segment>> segments, const std::string& name, int start) {
	std::unique_ptr<Molecule> molecule = std::make_unique<Molecule>(&input, lat, segments, name);
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
	bool ischarged=false;
	ForEachNode([&](const Node& node, int) {
		const auto& seg = *Seg[node.segment];
		if (!seg.state_name.empty()) {
			for (const Real valence : seg.state_valence) {
				if (valence != 0) ischarged = true;
			}
		} else if (seg.valence != 0) {
			ischarged = true;
		}
	});
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
		const int seg = SegmentIndices()[0];
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
	const auto monomers = SegmentIndices();
	int profile = 1;
	for (int seg : monomers) {
		OUTPUT["phi_" + Seg[seg]->name] = {{"profile", profile++}};
	}
}

std::span<Real> Molecule::GetPointer(int profile) {
NAMICS_DBG("GetPointer for Mol " + name << std::endl);
	const int M = lat->M;
	const size_t Msize = M;
	if (profile == 0) {
		lat->set_bounds(phitot.data());
		return phitot;
	}
	const auto monomers = SegmentIndices();
	if (profile > 0 && profile <= static_cast<int>(monomers.size())) {
		auto data = std::span<Real>(phi).subspan(static_cast<size_t>(profile - 1) * Msize, Msize);
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

void Molecule::propagate_backward(Real* G1, int &s, int block, int M, std::span<Real> ranked_phi) {
NAMICS_DBG("propagate_backward for Mol " + name << std::endl);

	int N= n_mon[block];
	const size_t Msize = M;
	for (int k=0; k<N; k++) {
		if (s<chainlength-1) {
			lat->propagate(Gg_b.data(),G1,(s+1)%2,s%2,M);
		} else {
			lat->Initiate(Gg_b.data()+(s%2)*M,G1);
		}

		lat->AddPhiS(phi.data()+molmon_nr[block]*M, Gg_f.data()+(s*M), Gg_b.data()+(s%2)*M);
		if (!ranked_phi.empty()) lat->AddPhiS(ranked_phi.data()+static_cast<size_t>(s)*Msize, Gg_f.data()+(s*M), Gg_b.data()+(s%2)*M);
		s--;
	}
}

bool Molecule::ComputePhi(){
NAMICS_DBG("ComputePhi for Molecule " + name << std::endl);
	return ComputePhiRanked({});
}


Real Molecule::fraction(int segnr){
NAMICS_DBG("fraction for Molecule " + name << std::endl); // Default for a monomer.
	int Nseg=0;
	ForEachNode([&](const Node& node, int) {
		if (segnr == node.segment) Nseg += node.repeat;
	});
	return 1.0*Nseg/chainlength;
}

Real* Molecule::ForwardBranch(int generation, int &s) {
NAMICS_DBG("ForwardBranch in Molecule " << std::endl);
	int b0 = first_b[generation];
	int bN = last_b[generation];
	std::vector<int> Br;
	std::vector<std::span<Real>> Gb;
	int M=lat->M;
	const size_t Msize = M;
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
					Gb.push_back(std::span<Real>(ForwardBranch(Gnr[k],s), Msize));
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

void Molecule::BackwardBranch(int generation, int &s, std::span<Real> ranked_phi){
NAMICS_DBG("BackwardBranch in Molecule " << std::endl);

	int b0 = first_b[generation];
	int bN = last_b[generation];
	std::vector<int> Br;
	std::vector<std::span<Real>> Gb;
	int M=lat->M;
	const size_t Msize = M;
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
				Gb.push_back(gg_f.subspan(static_cast<size_t>(last_s[Gnr[k]]) * Msize, Msize));
				ss=first_s[Gnr[k]];
				k-=(last_b[Gnr[k]]-first_b[Gnr[k]]+1) ;
			}
			Br.push_back(generation); ss--;
			Gb.push_back(gg_f.subspan(static_cast<size_t>(ss) * Msize, Msize));
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
					BackwardBranch(Br[i],s,ranked_phi);
				}
			}
			k++;
		} else {
			propagate_backward(Seg[mon_nr[k]]->G1.data(),s,k,M,ranked_phi);
		}

	}
}

bool Molecule::ComputePhiRanked(std::span<Real> ranked_phi) {
NAMICS_DBG("ComputePhi in Molecule " << std::endl);

	int M=lat->M;
	if (chainlength == 1) {
		std::copy_n(Seg[mon_nr[0]]->G1.begin(), M, phi.begin());
		GN = lat->WeightedSum(phi.data());
		for (int __i = 0; __i < M; ++__i) (phi)[__i] = (phi)[__i] * (Seg[mon_nr[0]]->G1)[__i];
		if (!ranked_phi.empty()) std::copy_n(phi.begin(), M, ranked_phi.begin());
		return true;
	}
	if (last_b.size() == 1) {
		int b0 = first_b[0];
		int bN = last_b[0];
		int s=0;
		Real* Glast=NULL;
		for (int b = b0; b<=bN ; ++b) Glast=propagate_forward(Seg[mon_nr[b]]->G1.data(),s,b,0,M);

		GN=lat->ComputeGN(Glast,M);

		s--;
		for (int b = bN ; b >= b0 ; --b) propagate_backward(Seg[mon_nr[b]]->G1.data(),s,b,M,ranked_phi);
		return true;
	}

	int generation=0;
	int s=0;
	Real* G=ForwardBranch(generation,s);
	GN=lat->ComputeGN(G,M);
	s--;
	BackwardBranch(generation,s,ranked_phi);

	return true;
}
