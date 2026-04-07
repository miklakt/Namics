#include "molecule.h"

#include <cctype>

namespace {

const std::vector<std::string>& MoleculeKeys() {
	static const std::vector<std::string> keys = {"freedom", "composition", "theta", "phibulk", "n", "B"};
	return keys;
}

using Topology = Molecule::Topology;
using Node = Molecule::Node;

struct MoleculeConfig {
	Molecule::OutputRequest output_request;
	std::string composition;
	std::string freedom;
	bool has_theta = false;
	bool has_n = false;
	bool has_phibulk = false;
	bool has_B = false;
	Real theta = 0;
	Real n = 0;
	Real phibulk = 0;
	Real B = 1;
};

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

bool ParseChain(const Molecule& mol, const std::string& s, size_t& pos, char terminator, Topology& chain, const std::unordered_map<std::string, int>& lookup);

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

bool ParseChain(const Molecule& mol, const std::string& s, size_t& pos, char terminator, Topology& chain, const std::unordered_map<std::string, int>& lookup) {
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
			++pos;
			Topology branch;
			if (!ParseChain(mol, s, pos, ']', branch, lookup)) {
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
			if (!ParseChain(mol, s, pos, ')', group, lookup)) {
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

bool ParseTopology(Molecule& mol, const std::string& composition) {
	mol.topology.clear();
	mol.segment_path.clear();
	mol.segment_types.clear();
	mol.chainlength = 0;

	std::unordered_map<std::string, int> lookup;
	for (size_t i = 0; i < mol.Seg.size(); ++i) lookup[mol.Seg[i]->name] = static_cast<int>(i);

	size_t pos = 0;
	if (!ParseChain(mol, composition, pos, '\0', mol.topology, lookup)) return false;
	SkipWhitespace(composition, pos);
	if (pos != composition.size() || mol.topology.empty()) {
		std::cout << "In composition of mol '" + mol.name + "' an invalid token was found: " << composition.substr(pos) << std::endl;
		return false;
	}
	return true;
}

bool BuildComposition(Molecule& mol, const std::string& composition){
NAMICS_DBG("BuildComposition for Mol " + mol.name << std::endl);
	if (!ParseTopology(mol, composition)) return false;
	std::vector<int> segment_type_index(mol.Seg.size(), -1);
	mol.segment_types.reserve(mol.Seg.size());
	const auto add_occurrence = [&](int segment, int parent) {
		if (segment_type_index[segment] < 0) {
			if (mol.Seg[segment]->freedom == "frozen") {
				std::cout << "In 'composition of mol " + mol.name + ", a segment was found with freedom 'frozen'. This is not permitted. " << std::endl;
				return -1;
			}
			segment_type_index[segment] = static_cast<int>(mol.segment_types.size());
			mol.segment_types.push_back(segment);
		}
		mol.segment_path.push_back({segment, segment_type_index[segment], {}});
		const int index = static_cast<int>(mol.segment_path.size()) - 1;
		if (parent >= 0) mol.segment_path[parent].children.push_back(index);
		return index;
	};
	const auto append_chain = [&](const auto& self, const Topology& chain, size_t pos, int parent) -> bool {
		if (pos >= chain.size()) return true;
		const Node& node = chain[pos];
		int tail = parent;
		for (int repeat = 0; repeat < node.repeat; ++repeat) {
			tail = add_occurrence(node.segment, tail);
			if (tail < 0) return false;
		}
		if (!self(self, chain, pos + 1, tail)) return false;
		for (const auto& branch : node.branches) {
			if (!self(self, branch, 0, tail)) return false;
		}
		return true;
	};
	if (!append_chain(append_chain, mol.topology, 0, -1)) return false;
	mol.chainlength = static_cast<int>(mol.segment_path.size());
	return mol.chainlength > 0;
}

bool ParseConfig(const Input& input, const std::string& name, int start, MoleculeConfig& config) {
	const auto& parameters = input.Parameters("mol", name, start);
	bool success = true;
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(MoleculeKeys(), it.key())) continue;
		success = false;
		std::cout << "mol property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : MoleculeKeys()) std::cout << item << std::endl;
	}
	if (!success) return false;

	try {
		config.composition = parameters.value("composition", std::string{});
		config.freedom = parameters.value("freedom", std::string{});
		config.has_theta = parameters.contains("theta");
		config.has_n = parameters.contains("n");
		config.has_phibulk = parameters.contains("phibulk");
		config.has_B = parameters.contains("B");
		if (config.has_theta) config.theta = parameters.at("theta").get<Real>();
		if (config.has_n) config.n = parameters.at("n").get<Real>();
		if (config.has_phibulk) config.phibulk = parameters.at("phibulk").get<Real>();
		if (config.has_B) config.B = parameters.at("B").get<Real>();
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in mol '" << name << "': " << error.what() << std::endl;
		return false;
	}
	return true;
}

bool HasMolOutputProperty(const Input& input, const std::string& name, int start, const std::string& property) {
	const auto& problem = input.Start(start);
	if (!problem.is_object()) return false;
	const auto json_it = problem.find("json");
	if (json_it == problem.end() || !json_it->is_object()) return false;
	const auto mol_it = json_it->find("mol");
	if (mol_it == json_it->end() || !mol_it->is_object()) return false;
	const auto has_property = [&](ParameterStore::const_iterator it) {
		if (it == mol_it->end()) return false;
		const auto& value = it.value();
		if (value.is_string()) return value.get<std::string>() == property;
		if (value.is_array()) {
			for (const auto& entry : value) {
				if (entry.is_string() && entry.get<std::string>() == property) return true;
			}
		}
		return false;
	};
	return has_property(mol_it->find(name)) || has_property(mol_it->find("*"));
}

void ParseOutputRequests(const Input& input, const Molecule& mol, int start, MoleculeConfig& config) {
	config.output_request.ranked_density = HasMolOutputProperty(input, mol.name, start, "phi_ranked");
	for (int seg : mol.SegmentTypes()) {
		if (HasMolOutputProperty(input, mol.name, start, "phi_" + mol.Seg[seg]->name)) {
			config.output_request.segment_density = true;
			return;
		}
	}
}

bool ConfigureFromConfig(Molecule& mol, const MoleculeConfig& config) {
NAMICS_DBG("Molecule:: ConfigureFromConfig for mol " << mol.name << std::endl);
	mol.composition = config.composition;
	mol.phibulk = 0;
	mol.n = 0;
	mol.theta = 0;
	mol.norm = 0;
	if (config.composition.empty()) {
		std::cout << "For mol '" + mol.name + "' the definition of 'composition' is required" << std::endl;
		return false;
	}
	if (!BuildComposition(mol, config.composition)) {
		std::cout << "For mol '" + mol.name + "' the composition is rejected. " << std::endl;
		return false;
	}
	const bool pinned = mol.IsPinned();
	if (config.freedom.empty()) {
		if (pinned) {
			std::cout <<"For mol " + mol.name + " the setting for 'freedom' was not set" << std::endl;
			return false;
		}
		std::cout <<"For mol " + mol.name + " the setting 'freedom' is expected: options: 'free' 'restricted' 'solvent' 'neutralizer' . Problem terminated " << std::endl;
		return false;
	}
	const bool allowed_freedom = config.freedom == "restricted" || (!pinned && (config.freedom == "free" || config.freedom == "solvent" || config.freedom == "neutralizer"));
	if (!allowed_freedom) {
		std::cout << "In mol " + mol.name + " the value for 'freedom' is not recognised " << std::endl;
		std::cout << "Select from: " << std::endl;
		if (!pinned) std::cout << "free ; solvent ; neutralizer ; ";
		std::cout << "restricted ; " << std::endl;
		return false;
	}
	mol.freedom = config.freedom;
	if (mol.freedom == "neutralizer" && !mol.IsCharged()) {
		std::cout << "Mol '" + mol.name + "' is not 'charged' and therefore this molecule can not be the neutralizer" << std::endl;
		return false;
	}
	if (mol.freedom == "free") {
		if (!config.has_phibulk) {
			std::cout <<"In mol " + mol.name + ", the setting 'freedom = free' should be combined with a value for 'phibulk'. "<<std::endl;
			return false;
		}
		mol.phibulk = config.phibulk;
		if (mol.phibulk < 0 || mol.phibulk > 1) {
			std::cout << "In mol " + mol.name + ", the value of 'phibulk' is out of range 0 .. 1." << std::endl;
			return false;
		}
	}
	mol.B = 1;
	if (!pinned && config.has_B) {
		mol.B = config.B;
		if (mol.B < 1e-9) {
			std::cout <<"for Mol" + mol.name + " mobility B should have a posititve value. Default value B=1 is chosen. " << std::endl;
			mol.B = 1;
		}
	}
	if (mol.freedom == "restricted") {
		if (!config.has_theta && !config.has_n) {
			std::cout <<"In mol " + mol.name + ", the setting 'freedom = restricted' should be combined with a value for 'theta' or 'n'; do not use both settings! "<<std::endl;
			return false;
		}
		if (config.has_theta && config.has_n) {
			std::cout <<"In mol " + mol.name + ", the setting 'freedom = restricted' does not allow both 'n' and 'theta' "<<std::endl;
			return false;
		}
		if (config.has_n) {
			mol.n = config.n;
			mol.theta = mol.n * mol.chainlength;
		}
		if (config.has_theta) {
			mol.theta = config.theta;
			mol.n = mol.theta / mol.chainlength;
		}
		if (mol.theta < 0 || (!pinned && mol.theta > mol.lat->volume)) {
			std::cout << "In mol " + mol.name + ", the value of 'n' or 'theta' is out of range." << std::endl;
			return false;
		}
	}
	return true;
}

} // namespace

Molecule::Molecule(Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_, std::string name_)
	: name(name_), Seg(Seg_), lat(Lat_) {
NAMICS_DBG("Constructor for Mol " + name << std::endl);
}

void Molecule:: AllocateMemory() {
	NAMICS_DBG("AllocateMemory in Mol " + name << std::endl);
	const int M = lat->M;
	phi.clear();
	phi_ranked.clear();
	phitot.assign(M, 0);
	q_forward.assign(static_cast<size_t>(M) * static_cast<size_t>(chainlength), 0);
	G_unity.assign(M, 0);
}

bool Molecule:: PrepareForCalculations(std::span<const Real> KSAM) {
NAMICS_DBG("PrepareForCalculations in Mol " + name << std::endl);
	std::copy(KSAM.begin(), KSAM.end(), G_unity.begin());
	std::fill(phitot.begin(), phitot.end(), 0);
	return true;
}

void Molecule::FinalizeOutputs() {
	const size_t M = static_cast<size_t>(lat->M);
	phi.clear();
	phi_ranked.clear();
	if (!output_request.any()) return;
	if (output_request.segment_density) phi.assign(M * SegmentTypes().size(), 0);
	if (output_request.ranked_density) phi_ranked.assign(M * static_cast<size_t>(chainlength), 0);
	AccumulateDensity(false, {}, {});
}

bool Molecule::IsPinned() {
NAMICS_DBG("IsPinned for Mol " + name << std::endl);
	bool pinned = false;
	ForEachNode([&](const Node& node) { pinned |= Seg[node.segment]->freedom == "pinned"; });
	return pinned;
}

Real Molecule::Charge() {
NAMICS_DBG("Molecule:: Charge" << std::endl);
	Real charge=0;
	ForEachOccurrence([&](const Node& node) {
		const auto& seg = *Seg[node.segment];
		if (seg.state_name.size() > 1) {
			for (size_t i = 0; i < seg.state_name.size(); ++i) charge += seg.state_alphabulk[i] * seg.state_valence[i];
		} else {
			charge += seg.valence;
		}
	});
	return charge/chainlength;
}

namespace molecule_factory {

std::unique_ptr<Molecule> CreateChecked(const Input& input, Lattice* lat, std::span<const std::unique_ptr<Segment>> segments, const std::string& name, int start) {
	MoleculeConfig config;
	if (!ParseConfig(input, name, start, config)) return nullptr;
	std::unique_ptr<Molecule> molecule = std::make_unique<Molecule>(lat, segments, name);
	if (!ConfigureFromConfig(*molecule, config)) return nullptr;
	ParseOutputRequests(input, *molecule, start, config);
	molecule->output_request = config.output_request;
	return molecule;
}

} // namespace molecule_factory

bool Molecule::IsCharged() {
NAMICS_DBG("IsCharged for Mol " + name << std::endl);
	bool ischarged=false;
	ForEachNode([&](const Node& node) {
		const auto& seg = *Seg[node.segment];
		if (!seg.state_name.empty()) {
			for (const Real valence : seg.state_valence) ischarged |= valence != 0;
		} else if (seg.valence != 0) {
			ischarged = true;
		}
	});
	return ischarged;
}

void Molecule::PushOutput() {
NAMICS_DBG("PushOutput for Mol " + name << std::endl);
	OUTPUT = nlohmann::ordered_json::object();
	OUTPUT["composition"] = composition;
	OUTPUT["freedom"] = freedom;
	if (freedom == "free") theta = lat->WeightedSum(phitot.data());
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
	if (chainlength==1) {
		const int seg = SegmentTypes()[0];
		if (Seg[seg]->ns >1) {
			if (mu_state.empty()) mu_state.assign(Seg[seg]->ns, Mu);
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
	const auto segment_types = SegmentTypes();
	int profile = 1;
	for (int seg : segment_types) {
		OUTPUT["phi_" + Seg[seg]->name] = {{"profile", profile++}};
	}
}

std::span<Real> Molecule::GetPointer(int profile) {
NAMICS_DBG("GetPointer for Mol " + name << std::endl);
	const int M = lat->M;
	const size_t Msize = M;
	if (profile == 0) return phitot;
	const auto segment_types = SegmentTypes();
	if (profile > 0 && profile <= static_cast<int>(segment_types.size())) {
		if (phi.empty()) return {};
		return std::span<Real>(phi).subspan(static_cast<size_t>(profile - 1) * Msize, Msize);
	}
	return {};
}

void Molecule::AddSegmentDensity(int segment_index, std::span<const Real> q_backward, bool include_segment_phi, std::span<Real> system_phitot, std::span<Real> mol_phitot) {
	const int M = lat->M;
	const auto& current = segment_path[segment_index];
	const int seg = current.segment;
	auto q = std::span<const Real>(q_forward).subspan(static_cast<size_t>(segment_index) * M, static_cast<size_t>(M));
	std::span<Real> mol_phi_block;
	if (!phi.empty()) {
		mol_phi_block = std::span<Real>(phi).subspan(static_cast<size_t>(current.segment_type_index) * M, static_cast<size_t>(M));
	}
	std::span<Real> ranked_phi_block;
	if (!phi_ranked.empty()) {
		ranked_phi_block = std::span<Real>(phi_ranked).subspan(static_cast<size_t>(segment_index) * M, static_cast<size_t>(M));
	}
	for (int i = 0; i < M; ++i) {
		Real value = q[i] * q_backward[i];
		if (norm > 0) value *= norm;
		if (include_segment_phi) Seg[seg]->phi[i] += value;
		if (!system_phitot.empty()) system_phitot[i] += value;
		if (!mol_phitot.empty()) mol_phitot[i] += value;
		if (!mol_phi_block.empty()) mol_phi_block[i] += value;
		if (!ranked_phi_block.empty()) ranked_phi_block[i] += value;
	}
}

void Molecule::PropagateForward() {
NAMICS_DBG("PropagateForward for Molecule " + name << std::endl);
	const int M = lat->M;
	std::vector<Real> work(2 * M);
	auto q_values = std::span<Real>(q_forward);
	for (int segment_index = chainlength - 1; segment_index >= 0; --segment_index) {
		const auto& current = segment_path[segment_index];
		auto q = q_values.subspan(static_cast<size_t>(segment_index) * M, static_cast<size_t>(M));
		lat->Initiate(q.data(), Seg[current.segment]->G1.data());
		for (int child : current.children) {
			auto q_child = q_values.subspan(static_cast<size_t>(child) * M, static_cast<size_t>(M));
			std::copy(q_child.begin(), q_child.end(), work.begin());
			lat->propagate(work.data(), G_unity.data(), 0, 1, M);
			for (int i = 0; i < M; ++i) q[i] *= work[M + i];
		}
	}
}

bool Molecule::ComputeGN(){
NAMICS_DBG("ComputeGN for Molecule " + name << std::endl);
	PropagateForward();
	const int M = lat->M;
	GN = lat->ComputeGN(q_forward.data(), M);
	return true;
}

void Molecule::AccumulateDensity(std::span<Real> system_phitot) {
	AccumulateDensity(true, system_phitot, std::span<Real>(phitot));
}


Real Molecule::fraction(int segnr){
NAMICS_DBG("fraction for Molecule " + name << std::endl); // Default for a monomer.
	int Nseg=0;
	ForEachOccurrence([&](const Node& node) {
		if (segnr == node.segment) ++Nseg;
	});
	return static_cast<Real>(Nseg)/chainlength;
}

void Molecule::PropagateBackward(int segment_index, std::span<const Real> q_backward, bool include_segment_phi, std::span<Real> system_phitot, std::span<Real> mol_phitot) {
	const int M = lat->M;
	std::vector<Real> current_q_backward(q_backward.begin(), q_backward.end());
	std::vector<Real> work(2 * M);
	while (segment_index >= 0) {
		AddSegmentDensity(segment_index, current_q_backward, include_segment_phi, system_phitot, mol_phitot);
		const auto& current = segment_path[segment_index];
		if (current.children.empty()) return;
		if (current.children.size() == 1) {
			const auto& g1 = Seg[current.segment]->G1;
			for (int i = 0; i < M; ++i) work[i] = current_q_backward[i] * g1[i];
			lat->propagate(work.data(), G_unity.data(), 0, 1, M);
			std::copy(work.begin() + M, work.end(), current_q_backward.begin());
			segment_index = current.children.front();
			continue;
		}

		const size_t child_count = current.children.size();
		std::vector<Real> shifted(child_count * static_cast<size_t>(M));
		const auto& g1 = Seg[current.segment]->G1;
		for (size_t child_index = 0; child_index < child_count; ++child_index) {
			auto q_child = std::span<const Real>(q_forward).subspan(static_cast<size_t>(current.children[child_index]) * M, static_cast<size_t>(M));
			std::copy(q_child.begin(), q_child.end(), work.begin());
			lat->propagate(work.data(), G_unity.data(), 0, 1, M);
			std::copy(work.begin() + M, work.end(), shifted.begin() + static_cast<std::ptrdiff_t>(child_index * static_cast<size_t>(M)));
		}

		const auto child_q_backward = [&](size_t excluded_child) {
			std::vector<Real> result(2 * M);
			for (int i = 0; i < M; ++i) result[i] = current_q_backward[i] * g1[i];
			for (size_t child_index = 0; child_index < child_count; ++child_index) {
				if (child_index == excluded_child) continue;
				const size_t base = child_index * static_cast<size_t>(M);
				for (int i = 0; i < M; ++i) result[i] *= shifted[base + static_cast<size_t>(i)];
			}
			lat->propagate(result.data(), G_unity.data(), 0, 1, M);
			return result;
		};

		for (size_t child_index = 1; child_index < child_count; ++child_index) {
			auto branch_q_backward = child_q_backward(child_index);
			PropagateBackward(current.children[child_index], std::span<const Real>(branch_q_backward).subspan(static_cast<size_t>(M), static_cast<size_t>(M)), include_segment_phi, system_phitot, mol_phitot);
		}
		auto q_continuation = child_q_backward(0);
		std::copy(q_continuation.begin() + M, q_continuation.end(), current_q_backward.begin());
		segment_index = current.children.front();
	}
}

bool Molecule::AccumulateDensity(bool include_segment_phi, std::span<Real> system_phitot, std::span<Real> mol_phitot) {
	if (segment_path.empty()) return true;
	std::vector<Real> q_backward(lat->M, 1.0);
	PropagateBackward(0, q_backward, include_segment_phi, system_phitot, mol_phitot);
	return true;
}
