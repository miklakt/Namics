#include "output.h"
#include "io_utils.h"

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

	json problem = {
		{"problem", start},
		{"name", base_name}
	};
	const int a = write_bounds ? 0 : lat->fjc;
	auto write_profile_json = [&](std::span<const Real> profile) {
		json out = json::array();
		switch (lat->gradients) {
			case 1:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; ++x) out.push_back(profile[x]);
				break;
			case 2:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; ++x) {
					json row = json::array();
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; ++y) row.push_back(profile[lat->P(x, y)]);
					out.push_back(std::move(row));
				}
				break;
			case 3:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; ++x) {
					json plane = json::array();
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; ++y) {
						json row = json::array();
						for (int z = a; z < lat->MZ + 2 * lat->fjc - a; ++z) row.push_back(profile[lat->P(x, y, z)]);
						plane.push_back(std::move(row));
					}
					out.push_back(std::move(plane));
				}
				break;
			default:
				break;
		}
		return out;
	};
	auto write_ranked_profile = [&](const Molecule& mol) {
		const int M = lat->M;
		const int ranks = mol.chainlength;
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
		return ranked;
	};
	bool has_profile_output = false;
	auto emit_output = [&](const std::string& key, const std::string& item_name, const std::string& item_prop, json value) {
		problem[key][item_name][item_prop] = std::move(value);
	};
	for (const auto& item : items) {
		const std::string key = item.value("key", "");
		const std::string item_name = item.value("name", "");
		const std::string item_prop = item.value("prop", "");
		std::string label = key;
		label.append(sep).append(item_name).append(sep).append(item_prop);

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
			emit_output(key, item_name, item_prop, nullptr);
			continue;
		}

		auto value_it = source->find(item_prop);
		if (value_it == source->end()) {
			std::cout << "Warning: unable to resolve json output quantity '" << label << "'" << std::endl;
			emit_output(key, item_name, item_prop, nullptr);
			continue;
		}

		std::span<const Real> profile;
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
				emit_output(key, item_name, item_prop, write_ranked_profile(*Mol[source_index]));
				continue;
			}
		}

		if (!profile.empty()) {
			emit_output(key, item_name, item_prop, write_profile_json(profile));
			has_profile_output = true;
			continue;
		}

		emit_output(key, item_name, item_prop, *value_it);
	};

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

	if (Sys->write_initial_guess) {
		json initial_guess;
		const std::span<const Real> values = std::span<const Real>(New->xx).first(static_cast<size_t>(New->iv));
		std::vector<std::string> monlist;
		std::vector<std::string> statelist;
		monlist.reserve(Sys->ItMonList.size());
		statelist.reserve(Sys->ItStateList.size());
		for (const int mon_index : Sys->ItMonList) monlist.push_back(Seg[mon_index]->name);
		for (const int state_index : Sys->ItStateList) statelist.push_back(Sta[state_index]->name);
		if (io::detail::WriteInitialGuessProfiles(initial_guess, values, monlist, statelist, Sys->charged, lat->M)) {
			problem["initial_guess"] = std::move(initial_guess);
		}
	}

	if (has_profile_output) {
		std::vector<Real> x_values;
		std::vector<Real> y_values;
		std::vector<Real> z_values;
		switch (lat->gradients) {
			case 1:
				x_values.reserve(static_cast<size_t>(lat->MX + 2 * lat->fjc - 2 * a));
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; ++x) x_values.push_back(coord_x(x));
				break;
			case 2:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; ++x) {
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; ++y) {
						x_values.push_back(coord_x(x));
						y_values.push_back(coord_y(y));
					}
				}
				break;
			case 3:
				for (int x = a; x < lat->MX + 2 * lat->fjc - a; ++x) {
					for (int y = a; y < lat->MY + 2 * lat->fjc - a; ++y) {
						for (int z = a; z < lat->MZ + 2 * lat->fjc - a; ++z) {
							x_values.push_back(coord_x(x));
							y_values.push_back(coord_y(y));
							z_values.push_back(coord_z(z));
						}
					}
				}
				break;
			default:
				break;
		}
		problem["x"] = std::move(x_values);
		if (lat->gradients >= 2) problem["y"] = std::move(y_values);
		if (lat->gradients >= 3) problem["z"] = std::move(z_values);
	}

	json metadata = {
		{"name", In->json_path}
	};

	const bool first_problem_of_run = (start == 1 && subl == 0);
	if (!json_writer.WriteProblem(filename, problem, metadata, append, first_problem_of_run)) {
		std::cout << "Failed to write json output file " << filename << std::endl;
	}
}
