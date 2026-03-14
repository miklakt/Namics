
#include "output.h"

#include <filesystem>
#include <limits>

namespace {

std::string JsonEscape(const std::string& value) {
	std::string out;
	out.reserve(value.size());
	for (const char c : value) {
		switch (c) {
			case '\\': out += "\\\\"; break;
			case '"': out += "\\\""; break;
			case '\n': out += "\\n"; break;
			case '\r': out += "\\r"; break;
			case '\t': out += "\\t"; break;
			default: out.push_back(c); break;
		}
	}
	return out;
}

template <typename T>
std::string JsonNumber(T value) {
	std::ostringstream out;
	out << std::setprecision(std::numeric_limits<T>::digits10 + 2) << value;
	return out.str();
}

bool IsJsonBoolString(const std::string& value) {
	return value == "true" || value == "false";
}

int GuessProfileSize(int mx, int my, int mz, int fjc) {
	if (my == 0) return (mx + 2 * fjc);
	if (mz == 0) return (mx + 2 * fjc) * (my + 2 * fjc);
	return (mx + 2 * fjc) * (my + 2 * fjc) * (mz + 2 * fjc);
}

std::string BuildInitialGuessObject(const std::string& method,
                                    int problem,
                                    int mx,
                                    int my,
                                    int mz,
                                    int fjc,
                                    bool charged,
                                    const std::vector<std::string>& monlist,
                                    const std::vector<std::string>& statelist,
                                    std::span<const Real> values) {
	const int m = GuessProfileSize(mx, my, mz, fjc);
	const int expected = static_cast<int>(monlist.size() + statelist.size() + (charged ? 1 : 0)) * m;
	if (static_cast<int>(values.size()) != expected) return "";

	std::ostringstream out;
	auto write_string_array = [&](const std::vector<std::string>& names) {
		out << "[";
		for (size_t i = 0; i < names.size(); ++i) {
			if (i > 0) out << ", ";
			out << "\"" << JsonEscape(names[i]) << "\"";
		}
		out << "]";
	};
	auto write_profile = [&](const std::string& key, std::span<const Real> data, bool trailing_comma) {
		out << "          \"" << JsonEscape(key) << "\": [";
		for (int i = 0; i < m; ++i) {
			if (i > 0) out << ", ";
			out << JsonNumber(data[i]);
		}
		out << "]";
		if (trailing_comma) out << ",";
		out << "\n";
	};

	out << "{\n";
	out << "        \"metadata\": {\n";
	out << "          \"problem\": " << problem << ",\n";
	out << "          \"method\": \"" << JsonEscape(method) << "\",\n";
	out << "          \"mx\": " << mx << ",\n";
	out << "          \"my\": " << my << ",\n";
	out << "          \"mz\": " << mz << ",\n";
	out << "          \"fjc\": " << fjc << ",\n";
	out << "          \"charged\": " << (charged ? "true" : "false") << "\n";
	out << "        },\n";
	out << "        \"monlist\": ";
	write_string_array(monlist);
	out << ",\n";
	out << "        \"statelist\": ";
	write_string_array(statelist);
	out << ",\n";
	out << "        \"profiles\": {\n";

	int offset = 0;
	const int total_profiles = static_cast<int>(monlist.size() + statelist.size() + (charged ? 1 : 0));
	int profile_index = 0;
	for (const std::string& mon_name : monlist) {
		++profile_index;
		write_profile("mon:" + mon_name, values.subspan(static_cast<size_t>(offset * m), static_cast<size_t>(m)), profile_index < total_profiles);
		++offset;
	}
	for (const std::string& state_name : statelist) {
		++profile_index;
		write_profile("state:" + state_name, values.subspan(static_cast<size_t>(offset * m), static_cast<size_t>(m)), profile_index < total_profiles);
		++offset;
	}
	if (charged) {
		++profile_index;
		write_profile("psi", values.subspan(static_cast<size_t>(offset * m), static_cast<size_t>(m)), profile_index < total_profiles);
	}
	out << "        }\n";
	out << "      }";
	return out.str();
}

} // namespace

Output::Output(const Input* In_,Lattice* Lat_,std::span<const std::unique_ptr<Segment>> Seg_,std::span<const std::unique_ptr<State>> Sta_, std::span<const std::unique_ptr<Reaction>> Rea_, std::span<const std::unique_ptr<Molecule>> Mol_,System* Sys_,Solve_scf* New_,std::string name_) {
NAMICS_DBG("constructor in Output "<< std::endl);	In=In_; Seg=Seg_; Sta=Sta_; Rea=Rea_; Mol=Mol_; Sys=Sys_; name=name_; New=New_;
	json_writer = io::json::SharedJsonWriter();
	lat=Lat_;
	KEYS.push_back("write_bounds");
	KEYS.push_back("append");
	KEYS.push_back("write");
	KEYS.push_back("header_separator");
	KEYS.push_back("filename");
}
Output::~Output() {
NAMICS_DBG("destructor in output " << std::endl);}
void Output::PutParameter(std::string new_param) {
NAMICS_DBG("PutParameter in Output " << std::endl); KEYS.push_back(new_param);
}

bool Output::Load() {
NAMICS_DBG("Load in output " << std::endl);	bool success=true;
	OUT_key.clear();
	OUT_name.clear();
	OUT_prop.clear();

	success = In->LoadItems(name, OUT_key, OUT_name, OUT_prop);

	if (success) {
		for (size_t i = 0; i < OUT_key.size(); ) {
			if (OUT_key[i]=="mol"){
				const std::string& wildcard = OUT_prop[i];
				const size_t star = wildcard.find('*');
				if (star == std::string::npos) {
					++i;
					continue;
				}
				if (star == 0 || wildcard.find('*', star + 1) != std::string::npos) {
					std::cout << "Mol wildcard output requires exactly one '*' and an explicit monomer prefix." << std::endl;
					return false;
				}
				const std::string prefix = wildcard.substr(0, star);
				const std::string suffix = wildcard.substr(star + 1);

				int molnr = -1;
				if (!ContainsValue(In->MolList, OUT_name[i], &molnr)) {
					std::cout << "Program error: output references unknown molecule '" << OUT_name[i] << "'." << std::endl;
					return false;
				}
				const int monlength = Mol[molnr]->MolMonList.size();
				for (int j=0; j<monlength; j++) {
					OUT_key.push_back(OUT_key[i]);
					OUT_name.push_back(OUT_name[i]);
					std::string s=prefix;
					s=s.append(Seg[Mol[molnr]->MolMonList[j]]->name);
					s=s.append(suffix);
					OUT_prop.push_back(s);
				}
				OUT_key.erase(OUT_key.begin()+i);
				OUT_name.erase(OUT_name.begin()+i);
				OUT_prop.erase(OUT_prop.begin()+i);
				continue;
			}
			++i;
		}

	}
	return success;
}

bool Output::CheckInput(int start_) {
NAMICS_DBG("CheckInput in output " << std::endl);	start=start_;
	bool success=true;
	success=In->CheckParameters("output",name,start, KEYS, PARAMETERS);
	if (success) {
		if (GetValue("append").size()>0) {
			append=ParseBool(GetValue("append"),append);
		} else {
			append=false;
		}

		write_bounds = ParseBool(GetValue("write_bounds"),false);
		write  = ParseBool(GetValue("write"),true);

		if (GetValue("header_separator").size()>0) {
			sep=GetValue("header_separator");
			if (sep=="classic") {
				sep = ":";
			} else {
				if (sep.length()>1) {
					std::cout <<"For output entry 'header_separator', expected to find the keyword 'classic' (meaning ':') or a single character." << std::endl;
					std::cout <<"header_separator set to the default value '_'" << std::endl;
					sep="_";
				}
			}
		} else {
			sep = "_";
		}

			if (success) {
				if (!Load()) {
					std::cout <<"Error in Load() in output" << std::endl;
				success=false;
			}
		}
	} else std::cout <<"Error in CheckParameters in output" << std::endl;
	return success;
}

std::string Output::GetValue(std::string parameter) {
NAMICS_DBG("GetValue in output " << std::endl); auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

std::span<int> Output::GetPointerInt(std::string key, std::string name, std::string prop) {
NAMICS_DBG("GetPointerInt in output " << std::endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int statelistlength=In->StateList.size();
	if (key=="output") return {};
	int listlength;
	int choice=0;
	int i,j;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if (key=="state") choice=4;
	if (key=="lat") choice=5;

	switch(choice) {
		case 1:
			listlength=Sys->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==Sys->strings[j]) return Sys->GetPointerInt(Sys->strings_value[j]);
				j++;
			}
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In->MolList[i]) {
					listlength= Mol[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Mol[i]->strings[j]) return Mol[i]->GetPointerInt(Mol[i]->strings_value[j]);
						j++;
					}
				}
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In->MonList[i]) {
					listlength= Seg[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Seg[i]->strings[j]) return Seg[i]->GetPointerInt(Seg[i]->strings_value[j]);
						j++;
					}
				}
				i++;
			}
			break;
		case 4:
			i=0;
			while (i<statelistlength){
				if (name==In->StateList[i]) {
					listlength= Sta[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Sta[i]->strings[j]) return Sta[i]->GetPointerInt(Sta[i]->strings_value[j]);
						j++;
					}
				}
				i++;
			}
			break;
		case 5:
			listlength= lat->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==lat->strings[j]) return lat->GetPointerInt(lat->strings_value[j]);
				j++;
			}
			break;
		default:
			std::cout << "Program error: in Output, GetPointerInt reaches default...." << std::endl;
	}
	return {};
}
std::span<Real> Output::GetPointer(std::string key, std::string name, std::string prop) {
NAMICS_DBG("GetPointer in output " << std::endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int statelistlength=In->StateList.size();
	if (key=="output") return {};
	int listlength;
	int choice=0;
	int i,j;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if (key=="state") choice=4;
	if (key=="lat") choice=5;

	switch(choice) {
		case 1:
			listlength=Sys->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==Sys->strings[j]) return Sys->GetPointer(Sys->strings_value[j]);
				j++;
			}
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In->MolList[i]) {
					listlength= Mol[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Mol[i]->strings[j]) return Mol[i]->GetPointer(Mol[i]->strings_value[j]);
						j++;
					}
				}
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In->MonList[i]) {
					listlength= Seg[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Seg[i]->strings[j]) return Seg[i]->GetPointer(Seg[i]->strings_value[j]);
						j++;
					}
				}
				i++;
			}
			break;
		case 4:
			i=0;
			while (i<statelistlength){
				if (name==In->StateList[i]) {
					if (prop=="phi") {
						return std::span<Real>(Seg[Sta[i]->mon_nr]->phi_state).subspan(static_cast<size_t>(Sta[i]->state_nr * lat->M), static_cast<size_t>(lat->M));
					}
					listlength= Sta[i]->strings.size();
					j=0;
					while (j<listlength) {
						if (prop==Sta[i]->strings[j]) return Sta[i]->GetPointer(Sta[i]->strings_value[j]);
						j++;
					}
				}
				i++;
			}
			break;
		case 5:
			listlength= lat->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==lat->strings[j]) return lat->GetPointer(lat->strings_value[j]);
				j++;
			}
			break;
		default:
			std::cout << "Program error: in Output, GetPointer reaches default...." << std::endl;
	}
	return {};
}
int Output::GetValue(std::string key, std::string name, std::string prop, int &int_result, Real &Real_result, std::string &string_result) {
NAMICS_DBG("GetValue (long) in output " << std::endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int statelistlength=In->StateList.size();
	int choice=0;
	int i;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if  (key=="state") choice=4;
	if  (key=="newton") choice=5;
	if  (key=="lat") choice=6;
	if  (key=="output") choice=7;
	switch(choice) {
		case 1:
			return Sys->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 2:
			i=0;
			while (i<mollistlength){
				if (name==In->MolList[i]) return Mol[i]->GetValue(prop,int_result,Real_result,string_result);
				i++;
			}
			break;
		case 3:
			i=0;
			while (i<monlistlength){
				if (name==In->MonList[i]) return Seg[i]->GetValue(prop,int_result,Real_result,string_result);
				i++;
			}
			break;
		case 4:
			i=0;
			while (i<statelistlength){
				if (name==In->StateList[i]) return Sta[i]->GetValue(prop,int_result,Real_result,string_result);
				i++;
			}
			break;
		case 5:
			return New->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 6:
			return lat->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 7:
			return GetValue(prop,name,int_result,Real_result,string_result);
			break;
		default:
			std::cout << "Program error: in Output, GetValue reaches default...." << std::endl;
	}
	return 0;
}

void Output::WriteOutput(int subl) {
NAMICS_DBG("WriteOutput in output " + name << std::endl);	lat->subl=subl;
	if (!write) return;
	std::string filename;
	std::string base_name;
	const std::string configured_filename = GetValue("filename");
	if (configured_filename.size() > 0) {
		std::filesystem::path configured_path(In->ResolvePath(configured_filename));
		if (configured_path.extension() != ".json") {
			configured_path.replace_extension(".json");
		}
		filename = configured_path.string();
		if (configured_path.has_stem()) {
			base_name = configured_path.stem().string();
		} else {
			base_name = configured_path.filename().string();
		}
	} else {
		std::filesystem::path out_path(In->name);
		out_path = out_path.filename();

		if (out_path.has_stem()) {
			base_name = out_path.stem().string();
		} else {
			base_name = out_path.filename().string();
		}
		if (base_name.empty()) {
			base_name = "output";
		}

		filename = base_name + ".json";
		filename = In->GetOutputPath() + filename;
	}

	std::vector<std::span<Real>> profile_pointer;
	std::vector<std::string> profile_header;
	std::vector<std::pair<std::string, std::string>> scalar_values;
	int length = OUT_key.size();
	for (int i=0; i<length; i++) {
		std::string label = OUT_key[i];
		label.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);

		std::vector<std::string> prop_sub;
		In->split(OUT_prop[i],'(',prop_sub);
		const bool indexed_scalar = (prop_sub.size() > 1 && prop_sub[0] != OUT_prop[i]);

		if (!indexed_scalar) {
			auto profile = GetPointer(OUT_key[i], OUT_name[i], OUT_prop[i]);
			if (!profile.empty()) {
				profile_pointer.push_back(profile);
				profile_header.push_back(label);
				continue;
			}
		}

		int int_result = 0;
		Real Real_result = 0;
		std::string string_result;
		const std::string value_key = prop_sub.size() > 0 ? prop_sub[0] : OUT_prop[i];
		const int result_nr = GetValue(OUT_key[i], OUT_name[i], value_key, int_result, Real_result, string_result);
		std::string literal = "null";
		if (result_nr == 1) {
			literal = JsonNumber(int_result);
		} else if (result_nr == 2) {
			literal = JsonNumber(Real_result);
		} else if (result_nr == 3) {
			if (indexed_scalar) {
				auto profile = GetPointer(OUT_key[i], OUT_name[i], value_key);
				if (!profile.empty()) {
					literal = JsonNumber(lat->GetValue(profile, prop_sub[1]));
				}
			} else if (IsJsonBoolString(string_result)) {
				literal = string_result;
			} else {
				literal = "\"" + JsonEscape(string_result) + "\"";
			}
		} else {
			std::cout << "Warning: unable to resolve json output quantity '" << label << "'" << std::endl;
		}
		scalar_values.push_back({label, literal});
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
			for (size_t i = 0; i < profile_pointer.size(); ++i) {
				column_values[c++].push_back(profile_pointer[i][index]);
			}
		};

		if (!profile_pointer.empty()) {
			switch (lat->gradients) {
				case 1:
					for (int x = a; x < lat->MX + 2 * lat->fjc - a; x++) {
						append_row(coord_x(x), 0, 0, x);
					}
					break;
				case 2:
					for (int x = a; x < lat->MX + 2 * lat->fjc - a; x++) {
						for (int y = a; y < lat->MY + 2 * lat->fjc - a; y++) {
							append_row(coord_x(x), coord_y(y), 0, lat->P(x, y));
						}
					}
					break;
				case 3:
					for (int x = a; x < lat->MX + 2 * lat->fjc - a; x++) {
						for (int y = a; y < lat->MY + 2 * lat->fjc - a; y++) {
							for (int z = a; z < lat->MZ + 2 * lat->fjc - a; z++) {
								append_row(coord_x(x), coord_y(y), coord_z(z), lat->P(x, y, z));
							}
						}
					}
					break;
				default:
					break;
			}
		}

		std::string initial_guess_object;
		if (Sys->write_initial_guess) {
			std::vector<std::string> guess_monlist;
			std::vector<std::string> guess_statelist;
			const int mon_length = Sys->ItMonList.size();
			const int state_length = Sys->ItStateList.size();
			guess_monlist.reserve(mon_length);
			guess_statelist.reserve(state_length);
			for (int i = 0; i < mon_length; ++i) guess_monlist.push_back(Seg[Sys->ItMonList[i]]->name);
			for (int i = 0; i < state_length; ++i) guess_statelist.push_back(Sta[Sys->ItStateList[i]]->name);

			initial_guess_object = BuildInitialGuessObject(New->SCF_method,
			                                               start,
			                                               lat->MX,
			                                               lat->MY,
			                                               lat->MZ,
			                                               lat->fjc,
			                                               Sys->charged,
			                                               guess_monlist,
			                                               guess_statelist,
			                                               std::span<const Real>(New->xx).first(static_cast<size_t>(New->iv)));
			if (initial_guess_object.empty()) {
				std::cout << "Warning: unable to serialize embedded initial guess for problem " << start << std::endl;
			}
		}

		std::ostringstream problem;
		problem << "    {\n";
		problem << "      \"problem\": " << start << ",\n";
		problem << "      \"name\": \"" << JsonEscape(base_name) << "\"";

		for (size_t i = 0; i < scalar_values.size(); ++i) {
			problem << ",\n";
			problem << "      \"" << JsonEscape(scalar_values[i].first) << "\": " << scalar_values[i].second;
		}

		for (size_t i = 0; i < column_names.size(); ++i) {
			problem << ",\n";
			problem << "      \"" << JsonEscape(column_names[i]) << "\": [";
			for (size_t j = 0; j < column_values[i].size(); ++j) {
				if (j > 0) problem << ", ";
				problem << JsonNumber(column_values[i][j]);
			}
			problem << "]";
		}
		if (!initial_guess_object.empty()) {
			problem << ",\n";
			problem << "      \"initial_guess\": " << initial_guess_object;
		}
		problem << "\n";
		problem << "    }";

		std::ostringstream metadata;
		metadata << "{\n";
		metadata << "    \"name\": \"" << JsonEscape(In->name) << "\"\n";
		metadata << "  }";

	const bool first_problem_of_run = (start == 1 && subl == 0);
	if (!json_writer->WriteProblem(filename, problem.str(), metadata.str(), append, first_problem_of_run)) {
		std::cout << "Failed to write json output file " << filename << std::endl;
	}
}

int Output::GetValue(std::string prop, std::string mod, int& int_result, Real& Real_result, std::string& string_result) {
	(void)mod;
  int i = 0;
  int_result=0;
  int length = ints.size();
  while (i < length) {
    if (prop == ints[i]) {
      int_result = ints_value[i];
      return 1;
    }
    i++;
  }
  i = 0;
  Real_result=0.0;
  length = Reals.size();
  while (i < length) {
    if (prop == Reals[i]) {
      Real_result = Reals_value[i];
      return 2;
    }
    i++;
  }
  i = 0;
  string_result="false";
  length = bools.size();
  while (i < length) {
    if (prop == bools[i]) {
      if (bools_value[i])
        string_result = "true";
      else
        string_result = "false";
      return 3;
    }
    i++;
  }
  i = 0;
  length = strings.size();
  while (i < length) {
    if (prop == strings[i]) {
      string_result = strings_value[i];
      return 3;
    }
    i++;
  }
  return 0;
}

void Output::push(std::string s, Real X) {
  Reals.push_back(s);
  Reals_value.push_back(X);
}
void Output::push(std::string s, int X) {
  ints.push_back(s);
  ints_value.push_back(X);
}
void Output::push(std::string s, bool X) {
  bools.push_back(s);
  bools_value.push_back(X);
}

void Output::push(std::string s, std::string X) {
  strings.push_back(s);
  strings_value.push_back(X);
}
