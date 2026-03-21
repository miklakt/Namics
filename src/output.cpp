#include "output.h"

#include <filesystem>
#include <nlohmann/json.hpp>

Output::Output(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_,std::span<const std::unique_ptr<State>> Sta_, std::span<const std::unique_ptr<Reaction>> Rea_, std::span<const std::unique_ptr<Molecule>> Mol_,System* Sys_,Solve_scf* New_,std::string name_) {
NAMICS_DBG("constructor in Output "<< std::endl);	In=In_; Seg=Seg_; Sta=Sta_; Rea=Rea_; Mol=Mol_; Sys=Sys_; name=name_; New=New_;
	lat=Lat_;
}

Output::~Output() {
NAMICS_DBG("destructor in output " << std::endl);}

bool Output::Load() {
NAMICS_DBG("Load in output " << std::endl);
	items = In->LoadItems(name);
	if (!items.is_array()) return false;
	nlohmann::ordered_json expanded = nlohmann::ordered_json::array();
	for (const auto& item : items) {
		if (item.value("key", "") != "mol") {
			expanded.push_back(item);
			continue;
		}
		const std::string wildcard = item.value("prop", "");
		const size_t star = wildcard.find('*');
		if (star == std::string::npos) {
			expanded.push_back(item);
			continue;
		}
		if (star == 0 || wildcard.find('*', star + 1) != std::string::npos) {
			std::cout << "Mol wildcard output requires exactly one '*' and an explicit monomer prefix." << std::endl;
			return false;
		}
		int molnr = -1;
		const std::string mol_name = item.value("name", "");
		if (!ContainsValue(In->MolList, mol_name, &molnr)) {
			std::cout << "Program error: output references unknown molecule '" << mol_name << "'." << std::endl;
			return false;
		}
		const std::string prefix = wildcard.substr(0, star);
		const std::string suffix = wildcard.substr(star + 1);
		for (size_t j=0; j<Mol[molnr]->MolMonList.size(); j++) {
			expanded.push_back({
				{"key", "mol"},
				{"name", mol_name},
				{"prop", prefix + Seg[Mol[molnr]->MolMonList[j]]->name + suffix}
			});
		}
	}
	items = std::move(expanded);
	return true;
}

bool Output::CheckInput(int start_) {
NAMICS_DBG("CheckInput in output " << std::endl);	start=start_;
	bool success=true;
	const auto& parameters = In->Parameters("output", name, start);
	static const std::vector<std::string> keys = {"write_bounds", "append", "write", "header_separator", "filename"};
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (ContainsValue(keys, it.key())) continue;
		success = false;
		std::cout << "output property '" << it.key() << "' is unknown. Select from: " << std::endl;
		for (const std::string& item : keys) std::cout << item << std::endl;
	}
	if (success) {
		try {
			append = parameters.value("append", false);
			write_bounds = parameters.value("write_bounds",false);
			write  = parameters.value("write",true);
			if (parameters.contains("filename")) (void)parameters.at("filename").get<std::string>();

			const std::string header_separator = parameters.value("header_separator", std::string{});
			if (!header_separator.empty()) {
				sep=header_separator;
				if (sep=="classic") {
					sep = ":";
				} else if (sep.length()>1) {
					std::cout <<"For output entry 'header_separator', expected to find the keyword 'classic' (meaning ':') or a single character." << std::endl;
					std::cout <<"header_separator set to the default value '_'" << std::endl;
					sep="_";
				}
			} else {
				sep = "_";
			}
		} catch (const nlohmann::json::exception& error) {
			std::cout << "Invalid json type in output '" << name << "': " << error.what() << std::endl;
			success = false;
		}

		if (success && !Load()) {
			std::cout <<"Error in Load() in output" << std::endl;
			success=false;
		}
	}
	return success;
}

void Output::WriteOutput(int subl) {
NAMICS_DBG("WriteOutput in output " + name << std::endl);	lat->subl=subl;
	if (!write) return;
	using json = nlohmann::ordered_json;
	std::string filename;
	std::string base_name;
	const auto& output_config = In->Parameters("output", name, start);
	const std::string configured_filename = output_config.value("filename", std::string{});
	if (configured_filename.size() > 0) {
		std::filesystem::path configured_path = In->ResolvePath(configured_filename);
		const std::string name = configured_path.filename().string();
		if (name.size() < 12 || name.compare(name.size() - 12, 12, ".output.json") != 0) {
			configured_path.replace_extension(".output.json");
		}
		filename = configured_path.string();
		base_name = configured_path.has_stem() ? configured_path.stem().string() : configured_path.filename().string();
	} else {
		std::filesystem::path out_path(In->json_path);
		out_path = out_path.filename();
		base_name = out_path.has_stem() ? out_path.stem().string() : out_path.filename().string();
		if (base_name.size() > 6 && base_name.compare(base_name.size() - 6, 6, ".input") == 0) base_name.resize(base_name.size() - 6);
		if (base_name.empty()) base_name = "output";
		filename = In->GetOutputPath() + base_name + ".output.json";
	}

	std::vector<std::span<Real>> profile_pointer;
	std::vector<std::string> profile_header;
	std::vector<std::pair<std::string, json>> scalar_values;
	auto write_ranked_profile = [&](const std::string& label, const Molecule& mol) {
		const int M = lat->M;
		const int ranks = mol.chainlength;
		const int a = write_bounds ? 0 : lat->fjc;
		json ranked = json::array();
		for (int r = 0; r < ranks; ++r) {
			json rank = json::array();
			for (int x = a; x < lat->MX + 2 * lat->fjc - a; ++x) {
				if (lat->gradients == 1) {
					rank.push_back(mol.phi_ranked[static_cast<size_t>(r) * M + x]);
					continue;
				}
				if (lat->gradients == 2) {
					json row = json::array();
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; ++y) row.push_back(mol.phi_ranked[static_cast<size_t>(r) * M + lat->P(x, y)]);
					rank.push_back(std::move(row));
					continue;
				}
				if (lat->gradients == 3) {
					json plane = json::array();
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; ++y) {
						json row = json::array();
						for (int z = a; z < lat->MZ + 2 * lat->fjc - a; ++z) row.push_back(mol.phi_ranked[static_cast<size_t>(r) * M + lat->P(x, y, z)]);
						plane.push_back(std::move(row));
					}
					rank.push_back(std::move(plane));
				}
			}
			ranked.push_back(std::move(rank));
		}
		scalar_values.push_back({label, std::move(ranked)});
	};
	json restart = {{"method", New->SCF_method},
	                {"mx", lat->MX},
	                {"my", lat->MY},
	                {"mz", lat->MZ},
	                {"fjc", lat->fjc},
	                {"charged", false},
	                {"monlist", json::array()},
	                {"statelist", json::array()},
	                {"profiles", json::object()}};
	for (const auto& item : items) {
		const std::string key = item.value("key", "");
		const std::string item_name = item.value("name", "");
		const std::string item_prop = item.value("prop", "");
		std::string label = key;
		label.append(sep).append(item_name).append(sep).append(item_prop);
		const std::string& value_key = item_prop;

		const ParameterStore* source = nullptr;
		int source_index = -1;
		if (key == "output") source = &output_config;
		if (key == "sys") source = &Sys->OUTPUT;
		if (key == "newton") source = &New->OUTPUT;
		if (key == "lat") source = &lat->OUTPUT;
		if (key == "mol" && ContainsValue(In->MolList, item_name, &source_index)) source = &Mol[source_index]->OUTPUT;
		if (key == "mon" && ContainsValue(In->MonList, item_name, &source_index)) source = &Seg[source_index]->OUTPUT;
		if (key == "state" && ContainsValue(In->StateList, item_name, &source_index)) source = &Sta[source_index]->OUTPUT;
		if (key == "reaction" && ContainsValue(In->ReactionList, item_name, &source_index)) source = &Rea[source_index]->OUTPUT;
		if (source == nullptr) {
			std::cout << "Warning: unable to resolve json output quantity '" << label << "'" << std::endl;
			scalar_values.push_back({label, nullptr});
			continue;
		}

		auto value_it = source->find(value_key);
		if (value_it == source->end()) {
			std::cout << "Warning: unable to resolve json output quantity '" << label << "'" << std::endl;
			scalar_values.push_back({label, nullptr});
			continue;
		}

		std::span<Real> profile;
		if (value_it->is_object() && value_it->contains("profile")) {
			const int profile_id = value_it->at("profile").get<int>();
			if (key == "sys") profile = Sys->GetPointer(profile_id);
			if (key == "lat") profile = lat->GetPointer(profile_id);
			if (key == "mol" && source_index >= 0) profile = Mol[source_index]->GetPointer(profile_id);
			if (key == "mon" && source_index >= 0) profile = Seg[source_index]->GetPointer(profile_id);
			if (key == "state" && source_index >= 0) profile = Sta[source_index]->GetPointer(profile_id);
		}
		if (value_it->is_object() && value_it->contains("ranked_profile")) {
			if (key == "mol" && source_index >= 0 && item_prop == "phi_ranked") {
				write_ranked_profile(label, *Mol[source_index]);
				continue;
			}
		}

		if (!profile.empty()) {
			if (key == "mon" && item_prop == "u") {
				restart["monlist"].push_back(item_name);
				restart["profiles"]["mon:" + item_name] = std::vector<Real>(profile.begin(), profile.end());
			}
			if (key == "state" && item_prop == "u") {
				restart["statelist"].push_back(item_name);
				restart["profiles"]["state:" + item_name] = std::vector<Real>(profile.begin(), profile.end());
			}
			if (key == "sys" && item_prop == "psi") {
				restart["charged"] = true;
				restart["profiles"]["psi"] = std::vector<Real>(profile.begin(), profile.end());
			}
			profile_pointer.push_back(profile);
			profile_header.push_back(label);
			continue;
		}

		scalar_values.push_back({label, *value_it});
	}

	const int a = write_bounds ? 0 : lat->fjc;
	const Real inv_fjc = static_cast<Real>(1.0) / static_cast<Real>(lat->fjc);
	const auto coord_x = [&](int x) -> Real {
		return lat->offset_first_layer * inv_fjc + static_cast<Real>(x - lat->fjc + 1) * inv_fjc - static_cast<Real>(0.5) * inv_fjc;
	};
	const auto coord_y = [&](int y) -> Real {
		return static_cast<Real>(y - lat->fjc + 1) * inv_fjc - static_cast<Real>(0.5) * inv_fjc;
	};
	const auto coord_z = [&](int z) -> Real {
		return static_cast<Real>(z - lat->fjc + 1) * inv_fjc - static_cast<Real>(0.5) * inv_fjc;
	};

	std::vector<std::string> column_names;
	std::vector<std::vector<Real>> column_values;
	if (!profile_pointer.empty()) {
		if (lat->gradients >= 1) {
			column_names.push_back("x");
			column_values.emplace_back();
		}
		if (lat->gradients >= 2) {
			column_names.push_back("y");
			column_values.emplace_back();
		}
		if (lat->gradients >= 3) {
			column_names.push_back("z");
			column_values.emplace_back();
		}
		for (size_t i = 0; i < profile_header.size(); ++i) {
			column_names.push_back(profile_header[i]);
			column_values.emplace_back();
		}
	}

	auto append_row = [&](Real x, Real y, Real z, int index) {
		size_t c = 0;
		if (lat->gradients >= 1) column_values[c++].push_back(x);
		if (lat->gradients >= 2) column_values[c++].push_back(y);
		if (lat->gradients >= 3) column_values[c++].push_back(z);
		for (size_t i = 0; i < profile_pointer.size(); ++i) column_values[c++].push_back(profile_pointer[i][index]);
	};

	if (!profile_pointer.empty()) {
		switch (lat->gradients) {
			case 1:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; x++) append_row(coord_x(x), 0, 0, x);
				break;
			case 2:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; x++) {
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; y++) append_row(coord_x(x), coord_y(y), 0, lat->P(x, y));
				}
				break;
			case 3:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; x++) {
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; y++) {
						for (int z = a; z < lat->MZ + 2 * lat->fjc - a; z++) append_row(coord_x(x), coord_y(y), coord_z(z), lat->P(x, y, z));
					}
				}
				break;
			default:
				break;
		}
	}

	json initial_guess;
	if (Sys->write_initial_guess) {
		std::vector<std::string> guess_monlist;
		std::vector<std::string> guess_statelist;
		const int mon_length = Sys->ItMonList.size();
		const int state_length = Sys->ItStateList.size();
		guess_monlist.reserve(mon_length);
		guess_statelist.reserve(state_length);
		for (int i = 0; i < mon_length; ++i) guess_monlist.push_back(Seg[Sys->ItMonList[i]]->name);
		for (int i = 0; i < state_length; ++i) guess_statelist.push_back(Sta[Sys->ItStateList[i]]->name);
		const std::span<const Real> values = std::span<const Real>(New->xx).first(static_cast<size_t>(New->iv));
		const int m = lat->MY == 0 ? (lat->MX + 2 * lat->fjc) :
		             lat->MZ == 0 ? (lat->MX + 2 * lat->fjc) * (lat->MY + 2 * lat->fjc) :
		                            (lat->MX + 2 * lat->fjc) * (lat->MY + 2 * lat->fjc) * (lat->MZ + 2 * lat->fjc);
		const int expected = static_cast<int>(guess_monlist.size() + guess_statelist.size() + (Sys->charged ? 1 : 0)) * m;
		if (static_cast<int>(values.size()) != expected) {
			std::cout << "Warning: unable to serialize embedded initial guess for problem " << start << std::endl;
		} else {
			initial_guess["method"] = New->SCF_method;
			initial_guess["mx"] = lat->MX;
			initial_guess["my"] = lat->MY;
			initial_guess["mz"] = lat->MZ;
			initial_guess["fjc"] = lat->fjc;
			initial_guess["charged"] = Sys->charged;
			initial_guess["monlist"] = guess_monlist;
			initial_guess["statelist"] = guess_statelist;
			initial_guess["profiles"] = json::object();
			int offset = 0;
			for (const std::string& mon_name : guess_monlist) {
				initial_guess["profiles"]["mon:" + mon_name] = std::vector<Real>(values.begin() + offset * m, values.begin() + (offset + 1) * m);
				++offset;
			}
			for (const std::string& state_name : guess_statelist) {
				initial_guess["profiles"]["state:" + state_name] = std::vector<Real>(values.begin() + offset * m, values.begin() + (offset + 1) * m);
				++offset;
			}
			if (Sys->charged) initial_guess["profiles"]["psi"] = std::vector<Real>(values.begin() + offset * m, values.begin() + (offset + 1) * m);
		}
	}

	json problem = {
		{"problem", start},
		{"name", base_name}
	};
	const bool wrote_initial_guess = !initial_guess.is_null();
	for (size_t i = 0; i < scalar_values.size(); ++i) problem[scalar_values[i].first] = std::move(scalar_values[i].second);
	for (size_t i = 0; i < column_names.size(); ++i) problem[column_names[i]] = column_values[i];
	if (wrote_initial_guess) problem["initial_guess"] = std::move(initial_guess);
	if (!wrote_initial_guess && !restart["profiles"].empty()) {
		problem["method"] = std::move(restart["method"]);
		problem["mx"] = std::move(restart["mx"]);
		problem["my"] = std::move(restart["my"]);
		problem["mz"] = std::move(restart["mz"]);
		problem["fjc"] = std::move(restart["fjc"]);
		problem["charged"] = std::move(restart["charged"]);
		problem["monlist"] = std::move(restart["monlist"]);
		problem["statelist"] = std::move(restart["statelist"]);
		problem["profiles"] = std::move(restart["profiles"]);
	}

	json metadata = {
		{"name", In->json_path}
	};

	const bool first_problem_of_run = (start == 1 && subl == 0);
	if (!json_writer.WriteProblem(filename, problem, metadata, append, first_problem_of_run)) {
		std::cout << "Failed to write json output file " << filename << std::endl;
	}
}
