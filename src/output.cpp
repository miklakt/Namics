
#include "output.h"

#include <cctype>
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

} // namespace

Output::Output(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_,vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_,System* Sys_,Solve_scf* New_,string name_) {
NAMICS_DBG("constructor in Output "<< endl);	In=In_; Seg=Seg_; Sta=Sta_; Rea=Rea_; Mol=Mol_; Sys=Sys_; name=name_; New=New_;
	json_writer = io::json::SharedJsonWriter();
	lat=Lat_;
	KEYS.push_back("write_bounds");
	KEYS.push_back("append");
	KEYS.push_back("use_output_folder");
	KEYS.push_back("write");
	KEYS.push_back("header_separator");
	KEYS.push_back("filename");
	use_output_folder = true; // LINUX ONLY, when you remove this, add it as a default to its CheckInputs part.
}
Output::~Output() {
NAMICS_DBG("destructor in output " << endl);}
void Output::PutParameter(string new_param) {
NAMICS_DBG("PutParameter in Output " << endl); KEYS.push_back(new_param);
}

bool Output::Load() {
NAMICS_DBG("Load in output " << endl);	bool success=true;
	int molnr=0;
	OUT_key.clear();
	OUT_name.clear();
	OUT_prop.clear();

	success = In->LoadItems(name, OUT_key, OUT_name, OUT_prop);

	if (success) {
		int length=OUT_key.size();

		for (int i=0; i<length; i++) {
			if (OUT_key[i]=="mol"){
				vector<string> sub;
				In->split(OUT_prop[i],'*',sub);


				if (!(sub[0]==OUT_prop[i])){
					bool wildmon;
					if (sub[0] =="") wildmon=false; else wildmon=true;

					int k=0; int mollength=In->MolList.size();
					while (k<mollength) {
						if (In->MolList[k]==OUT_name[i]) molnr=k;
						k++;
					}
					if (wildmon) {
						int monlength=Mol[molnr]->MolMonList.size();
						for (int j=0; j<monlength; j++) {
							OUT_key.push_back(OUT_key[i]);
							OUT_name.push_back(OUT_name[i]);
							string s=sub[0];
							s=s.append(Seg[Mol[molnr]->MolMonList[j]]->name);
							s=s.append(sub[1]);
							OUT_prop.push_back(s);
						}
					} else {
						cout << "Alias-based wildcard output is not supported in this minimal build." << endl;
						return false;
					}
					OUT_key.erase(OUT_key.begin()+i);
					OUT_name.erase(OUT_name.begin()+i);
					OUT_prop.erase(OUT_prop.begin()+i);

				}

			}
		}

	}
	return success;
}

bool Output::CheckInput(int start_) {
NAMICS_DBG("CheckInput in output " << endl);	start=start_;
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
					cout <<"For output entry 'header_separator', expected to find the keyword 'classic' (meaning ':') or a single character." << endl;
					cout <<"header_separator set to the default value '_'" << endl;
					sep="_";
				}
			}
		} else {
			sep = "_";
		}

		if (GetValue("use_output_folder").size()>0) {
			use_output_folder = ParseBool(GetValue("use_output_folder"),use_output_folder);
		} // default is set in the constructor

		if (success) {
			if (!Load()) {
				cout <<"Error in Load() in output" << endl;
				success=false;
			}
		}
	} else cout <<"Error in CheckParameters in output" << endl;
	return success;
}

string Output::GetValue(string parameter) {
NAMICS_DBG("GetValue in output " << endl); auto it = PARAMETERS.find(parameter);
	if (it != PARAMETERS.end()) return it->second;
	return "";
}

int* Output::GetPointerInt(string key, string name, string prop, int &Size) {
NAMICS_DBG("GetPointerInt in output " << endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int listlength;
	int choice;
	int i,j;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if (key=="lat") choice=4;
	if (key=="output") choice = 5;

	switch(choice) {
		case 1:
			listlength=Sys->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==Sys->strings[j]) return Sys->GetPointerInt(Sys->strings_value[j],Size);
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
						if (prop==Mol[i]->strings[j]) return Mol[i]->GetPointerInt(Mol[i]->strings_value[j],Size);
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
						if (prop==Seg[i]->strings[j]) return Seg[i]->GetPointerInt(Seg[i]->strings_value[j],Size);
						j++;
					}
				}
				i++;
			}
			break;
		case 4:
			listlength= lat->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==lat->strings[j]) return lat->GetPointerInt(lat->strings_value[j],Size);
				j++;
			}
			break;
		case 5:
			listlength=PointerVectorInt.size();
			j=0;
			while (j<listlength) {
				if (prop==strings[j]){ Size=SizeVectorInt[j];  return PointerVectorInt[j];}
				j++;
			}
			break;
		default:
			cout << "Program error: in Output, GetPointerInt reaches default...." << endl;
	}
	return NULL;
}
Real* Output::GetPointer(string key, string name, string prop, int &Size) {
NAMICS_DBG("GetPointer in output " << endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int listlength;
	int choice;
	int i,j;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if (key=="lat") choice=4;
	if (key=="output") choice = 5;

	switch(choice) {
		case 1:
			listlength=Sys->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==Sys->strings[j]) return Sys->GetPointer(Sys->strings_value[j],Size);
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
						if (prop==Mol[i]->strings[j]) return Mol[i]->GetPointer(Mol[i]->strings_value[j],Size);
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
						if (prop==Seg[i]->strings[j]) return Seg[i]->GetPointer(Seg[i]->strings_value[j],Size);
						j++;
					}
				}
				i++;
			}
			break;
		case 4:
			listlength= lat->strings.size();
			j=0;
			while (j<listlength) {
				if (prop==lat->strings[j]) return lat->GetPointer(lat->strings_value[j],Size);
				j++;
			}
			break;
		case 5:
			listlength=PointerVectorReal.size();
			j=0;
			while (j<listlength) {
				if (prop==strings[j]){ Size=SizeVectorReal[j];  return PointerVectorReal[j];}
				j++;
			}
			break;
		default:
			cout << "Program error: in Output, GetPointer reaches default...." << endl;
	}
	return NULL;
}
int Output::GetValue(string key, string name, string prop, int &int_result, Real &Real_result, string &string_result) {
NAMICS_DBG("GetValue (long) in output " << endl); int monlistlength=In->MonList.size();
	int mollistlength=In->MolList.size();
	int choice=0;
	int i;
	if  (key=="sys") choice=1;
	if  (key=="mol") choice=2;
	if  (key=="mon") choice=3;
	if  (key=="newton") choice=4;
	if  (key=="lat") choice=5;
	if  (key=="output") choice=6;
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
			return New->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 5:
			return lat->GetValue(prop,int_result,Real_result,string_result);
			break;
		case 6:
			return GetValue(prop,name,int_result,Real_result,string_result);
			break;
		default:
			cout << "Program error: in Output, GetValue reaches default...." << endl;
	}
	return 0;
}

void Output::WriteOutput(int subl) {
NAMICS_DBG("WriteOutput in output " + name << endl);	lat->subl=subl;
	if (!write) return;
	int Size=0;
	string filename;
	const string configured_filename = GetValue("filename");
	const string infilename = configured_filename.size() > 0 ? configured_filename : In->name;
	std::filesystem::path out_path(infilename);
	if (use_output_folder) {
		out_path = out_path.filename();
	}

	string base_name;
	if (out_path.has_stem()) {
		base_name = out_path.stem().string();
	} else {
		base_name = out_path.filename().string();
	}
	if (base_name.empty()) {
		base_name = "output";
	}

	filename = base_name + ".json";
	filename = In->output_info.getOutputPath() + filename;

	vector<Real*> profile_pointer;
	vector<string> profile_header;
	vector<pair<string, string>> scalar_values;
	int length = OUT_key.size();
	for (int i=0; i<length; i++) {
		string label = OUT_key[i];
		label.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);

		vector<string> prop_sub;
		In->split(OUT_prop[i],'(',prop_sub);
		const bool indexed_scalar = (prop_sub.size() > 1 && prop_sub[0] != OUT_prop[i]);

		if (!indexed_scalar) {
			Real* profile = GetPointer(OUT_key[i], OUT_name[i], OUT_prop[i], Size);
			if (profile != NULL) {
				profile_pointer.push_back(profile);
				profile_header.push_back(label);
				continue;
			}
		}

		int int_result = 0;
		Real Real_result = 0;
		string string_result;
		const string value_key = prop_sub.size() > 0 ? prop_sub[0] : OUT_prop[i];
		const int result_nr = GetValue(OUT_key[i], OUT_name[i], value_key, int_result, Real_result, string_result);
		string literal = "null";
		if (result_nr == 1) {
			literal = JsonNumber(int_result);
		} else if (result_nr == 2) {
			literal = JsonNumber(Real_result);
		} else if (result_nr == 3) {
			if (indexed_scalar) {
				Real* profile = GetPointer(OUT_key[i], OUT_name[i], value_key, Size);
				if (profile != NULL) {
					literal = JsonNumber(lat->GetValue(profile, prop_sub[1]));
				}
			} else if (IsJsonBoolString(string_result)) {
				literal = string_result;
			} else {
				literal = "\"" + JsonEscape(string_result) + "\"";
			}
		} else {
			cout << "Warning: unable to resolve json output quantity '" << label << "'" << endl;
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

		vector<string> column_names;
		vector<vector<Real>> column_values;
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

		std::ostringstream problem;
		problem << "    {\n";
		problem << "      \"problem\": " << start << ",\n";
		// TODO(variate): reintroduce subloop metadata when the variate module is restored.
		problem << "      \"name\": \"" << JsonEscape(name) << "\"";

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
		problem << "\n";
		problem << "    }";

		std::ostringstream metadata;
		metadata << "{\n";
		metadata << "    \"name\": \"" << JsonEscape(In->name) << "\"\n";
		metadata << "  }";

	const bool first_problem_of_run = (start == 1 && subl == 0);
	if (!json_writer->WriteProblem(filename, problem.str(), metadata.str(), append, first_problem_of_run)) {
		cout << "Failed to write json output file " << filename << endl;
	}
}

int Output::GetValue(string prop, string mod, int& int_result, Real& Real_result, string& string_result) {
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

void Output::push(string s, Real X) {
  Reals.push_back(s);
  Reals_value.push_back(X);
}
void Output::push(string s, int X) {
  ints.push_back(s);
  ints_value.push_back(X);
}
void Output::push(string s, bool X) {
  bools.push_back(s);
  bools_value.push_back(X);
}

void Output::push(string s, string X) {
  strings.push_back(s);
  strings_value.push_back(X);
}
