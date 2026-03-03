
#include "output.h"

#include <cctype>
#include <limits>
#include <unordered_set>

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

Output::Output(const Input* In_,Lattice* Lat_,vector<Segment*> Seg_,vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_,System* Sys_,Solve_scf* New_,string name_,int outnr,int N_out) {
NAMICS_DBG("constructor in Output "<< endl);	In=In_; Lat = Lat_; Seg=Seg_; Sta=Sta_; Rea=Rea_; Mol=Mol_; Sys=Sys_; name=name_; n_output=N_out; output_nr=outnr;  New=New_;
	writer = io::legacy::SharedWriter();
	json_writer = io::json::SharedJsonWriter();
	//KEYS.push_back("write_output");
	lat=Lat;
	KEYS.push_back("write_bounds");
	KEYS.push_back("append");
	KEYS.push_back("use_output_folder");
	KEYS.push_back("write");
	KEYS.push_back("clear");
	KEYS.push_back("DOS");
	KEYS.push_back("header_separator");
	KEYS.push_back("filename");
	input_error=false;
	bin_folder = "bin"; // folder in Namics where the binary is located
	use_output_folder = true; // LINUX ONLY, when you remove this, add it as a default to its CheckInputs part.
	n_starts = In->GetNumStarts();
	first=0;

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

	if (name == "json") {
		std::unordered_set<std::string> seen_entries;
		auto append_items = [&](const std::vector<std::string>& keys,
		                        const std::vector<std::string>& names,
		                        const std::vector<std::string>& props) {
			for (size_t i = 0; i < keys.size(); ++i) {
				const std::string signature = keys[i] + "\n" + names[i] + "\n" + props[i];
				if (seen_entries.insert(signature).second) {
					OUT_key.push_back(keys[i]);
					OUT_name.push_back(names[i]);
					OUT_prop.push_back(props[i]);
				}
			}
		};

		for (const std::string& template_name : {std::string("json"), std::string("kal"), std::string("pro")}) {
			std::vector<std::string> key;
			std::vector<std::string> out_name;
			std::vector<std::string> prop;
			if (!In->LoadItems(template_name, key, out_name, prop)) {
				success = false;
				break;
			}
			append_items(key, out_name, prop);
		}
	} else {
		success= In->LoadItems(name, OUT_key, OUT_name, OUT_prop);
	}

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
	if (name != "kal" && name != "pro" && name != "json") {
		write = false;
		cout << "Output type '" << name << "' is disabled. Only 'kal', 'pro' and 'json' are supported." << endl;
		return true;
	}

	bool success=true;
	success=In->CheckParameters("output",name,start, KEYS, PARAMETERS);
	if (success) {
		DOS=false;
		if (GetValue("DOS").size()>0) {
			DOS=ParseBool(GetValue("DOS"),DOS);
		}

			if (GetValue("append").size()>0) {
				append=ParseBool(GetValue("append"),append);

				if (name=="pro") {
						if (append) cout << "Warning: for output of type 'pro', the append is set to 'false'." << endl;
				}
				if (first==0) first=start;
			} else {
				if (name=="kal") append=false;
				if (name=="pro") append=false;
				if (name=="json") append=false;
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
		write_option="no_error";
		if (GetValue("write_output").size()>0) {
			vector<string> option_list;
			option_list.push_back("always");
			option_list.push_back("no_error");
			option_list.push_back("never");
			if (!ParseString(GetValue("write_output"),write_option,option_list,"In output: 'write_output' not recognised. Use 'always', 'never', or 'no_error'. The last value is default.")){
				cout <<"continue with write_output : no_error" << endl;
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
	if (name != "kal" && name != "pro" && name != "json") {
		return;
	}
	int Size=0;
	string filename;
	vector<string> sub;

	string infilename = In->name;
	if (GetValue("filename").size()>0) infilename=GetValue("filename");
	In->split(infilename,'.',sub);
	if (sub.size() == 0) sub.push_back(infilename);
	string key;

	if (use_output_folder == true) {

		int occurrences = 0;
		string::size_type path_pos = 0;
		string slash = "/";

		// Check if we have any slashes in the filename (are we in inputs, or higher up?)
		while ((path_pos = sub[0].find(slash, path_pos)) != string::npos) {
    	++occurrences;
    	path_pos += slash.length();
		}

		// If we're not in the inputs folder, discard the path and take only the filename
		if  (occurrences != 0) {
			size_t found = sub[0].find_last_of("/\\");
			sub[0] = sub[0].substr(found+1);
		}
	}

    string numc = to_string(subl);
    string numcc = to_string(start);

	if (name=="json") {
		filename = sub[0].append(".json");
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
		return;
	}

    if (name=="kal") filename=sub[0].append(".").append(name); else {
		if (n_starts==1 && subl < 1) filename=sub[0].append(".").append(name);
		if (n_starts==1 && subl >0) filename = sub[0].append("_").append(numc).append(".").append(name);
		if (n_starts>1  && subl < 1) filename = sub[0].append("_").append(numcc).append(".").append(name);
		if (n_starts>1 && subl >0)  filename=sub[0].append("_").append(numcc).append("_").append(numc).append(".").append(name);
	}

	filename = In->output_info.getOutputPath() + filename;
	if (name=="pro") {
		vector<Real*> pointer;
		FILE *fp;
		fp=writer->OpenRaw(filename.c_str(),"w");
		int length=OUT_key.size();
		switch(lat->gradients) {
			case 1:
				writer->Writef(fp,"x\t");
				break;
			case 2:
				writer->Writef(fp,"x\ty\t");
				break;
			case 3:
				writer->Writef(fp,"x\ty\tz\t");
				break;
			default:
				break;
		}
		for (int i=0; i<length; i++) {
			string s=OUT_key[i];
			s.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);
			Real*  X = GetPointer(OUT_key[i],OUT_name[i],OUT_prop[i],Size);
			if (X!=NULL) {
				pointer.push_back(X);
				key = OUT_key[i];
				s = key.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);
				if (i<length-1) writer->Writef(fp,"%s\t",s.c_str()); else writer->Writef(fp,"%s",s.c_str());
			} else {cout << " Error for 'pro' output. It is only possible to output quantities known to be a 'profile'. That is why output quantity " + s + " is rejected. " << endl;}
		}
		if (DOS) writer->Writef(fp,"\r\n"); else writer->Writef(fp,"\n");
		Lat-> PutProfiles(fp,pointer,write_bounds,DOS);

		writer->Close(fp);
	}

	if (name=="kal") {
		if (start>first || subl>0) append=true;
		FILE *fp;
		if (!(writer->Exists(filename) && append)) {
			fp=writer->OpenRaw(filename.c_str(),"w");
			int length = OUT_key.size();
			for (int i=0; i<length; i++) {
				key=OUT_key[i];
				string s=key.append(sep).append(OUT_name[i]).append(sep).append(OUT_prop[i]);
				if (i<length-1) writer->Writef(fp,"%s\t",s.c_str()); else writer->Writef(fp,"%s",s.c_str());
			}
			if (DOS) writer->Writef(fp,"\r\n"); else writer->Writef(fp,"\n");
		} else fp=writer->OpenRaw(filename.c_str(),"a");

		if (fp == NULL) {
			cerr << "Error trying to open " << filename.c_str() << endl;
			perror("Error");
		}


		int length = OUT_key.size();
		for (int i=0; i<length; i++) {
			int int_result=0;
			int result_nr=0;
			Real Real_result=0;
			string string_result;
			vector<string> sub;
			In-> split(OUT_prop[i],'(',sub);
			result_nr= GetValue(OUT_key[i],OUT_name[i],sub[0],int_result,Real_result,string_result);
			if (result_nr==0) {if (i<length-1) writer->Writef(fp,"NiN\t"); else writer->Writef(fp,"NiN");}
			if (result_nr==1) {if (i<length-1) writer->Writef(fp,"%i\t",int_result); else writer->Writef(fp,"%i",int_result);}
#ifdef LongReal
			if (result_nr==2) {if (i<length-1) writer->Writef(fp,"%.16Le\t",Real_result); else  writer->Writef(fp,"%.16Le",Real_result);}
#else
			if (result_nr==2) {if (i<length-1) writer->Writef(fp,"%.16e\t",Real_result); else  writer->Writef(fp,"%.16e",Real_result);}
#endif
			if (result_nr==3) {
				if (sub[0]==OUT_prop[i]) {
					if (i<length-1) writer->Writef(fp,"%s\t",string_result.c_str()); else writer->Writef(fp,"%s",string_result.c_str());
				} else {
					Real* X=GetPointer(OUT_key[i],OUT_name[i],sub[0],Size);
#ifdef LongReal
					if (i<length-1) writer->Writef(fp,"%.16Le\t",lat->GetValue(X,sub[1])); else writer->Writef(fp,"%.16Le",lat->GetValue(X,sub[1]));
#else
					if (i<length-1) writer->Writef(fp,"%.16e\t",lat->GetValue(X,sub[1])); else writer->Writef(fp,"%.16e",lat->GetValue(X,sub[1]));
#endif
				}
			}
		}
		if (DOS) writer->Writef(fp,"\r\n"); else writer->Writef(fp,"\n");
		writer->Close(fp);
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
