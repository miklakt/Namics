#include "input.h"
#include <cctype>
#include <filesystem>

namespace {
constexpr const char* OUTPUT_INFO_KEY = "out_info";
constexpr const char* DEFAULT_OUTPUT_PATH = "./output/";

void NormalizeLine(string& line, bool strip_tabs = true) {
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
	if (strip_tabs) {
		line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
	}
}

bool StartsWith(const string& value, const string& prefix) {
	return value.rfind(prefix, 0) == 0;
}

bool EndsWith(const string& value, const string& suffix) {
	return value.size() >= suffix.size() &&
		value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool IsEmptyOrComment(const string& line) {
	return line.empty() || (line.size() >= 2 && StartsWith(line, "//"));
}

string NormalizeOutputPath(string path) {
	if (path.empty()) return DEFAULT_OUTPUT_PATH;
	if (path.back() != '/') {
		path.push_back('/');
	}
	return path;
}

std::filesystem::path InputDirectoryFromName(const string& input_name) {
	std::filesystem::path input_path(input_name);
	const std::filesystem::path parent = input_path.parent_path();
	return parent.empty() ? std::filesystem::path(".") : parent;
}

bool HasBalancedBrackets(
	const string& expression,
	char open_bracket,
	char close_bracket,
	vector<int>& open_positions,
	vector<int>& close_positions) {
	vector<char> stack;
	for (size_t i = 0; i < expression.size(); ++i) {
		if (expression[i] == open_bracket) {
			stack.push_back(expression[i]);
			open_positions.push_back(static_cast<int>(i));
		} else if (expression[i] == close_bracket) {
			close_positions.push_back(static_cast<int>(i));
			if (stack.empty()) return false;
			stack.pop_back();
		}
	}
	return stack.empty();
}
} // namespace

Input::Input(const string& name_) {
	name=name_;
	output_path = DEFAULT_OUTPUT_PATH;
	KEYS = {
		"start",
		"sys",
		"mol",
		"mon",
		"lat",
		"newton",
		"output",
		OUTPUT_INFO_KEY,
		"state",
		"reaction"
	};

	in_file.open(name.c_str()); Input_error=false;

	if (in_file.is_open()) {
		int line_nr=0;
		std:: string In_line;
		std:: string last;
		while (in_file) {
			line_nr++;
			std::getline(in_file,In_line);
			NormalizeLine(In_line);
			bool add = !IsEmptyOrComment(In_line);
			if (In_line.length()>7 && StartsWith(In_line, "include")) {
				add = false;
				string filename_inc;
				const string include_spec = In_line.substr(8);
				if (EndsWith(include_spec, "::")) filename_inc = include_spec.substr(0, include_spec.size() - 2);
				else filename_inc = include_spec;
				filename_inc = ResolvePath(filename_inc);
				inc_file.open(filename_inc);
				if (inc_file.is_open()) {
					int line_nr_inc=0;
					while (inc_file) {
						line_nr_inc++;
						std::getline(inc_file,In_line);
						NormalizeLine(In_line);
						if (!IsEmptyOrComment(In_line)) {
							elems.push_back(std::to_string(line_nr).append("-").append(std::to_string(line_nr_inc)).append(":").append(In_line));
						}
					}

					inc_file.close();
				} else {
					cout <<"'include : " << filename_inc <<"' not working because file is not found " << endl; Input_error=true;
				}

			}
			if (add) {
				elems.push_back(std::to_string(line_nr).append(":").append(In_line));
				last=In_line;
			}
		}

		if (last.substr(0,5) != "start") {
		   elems.push_back(std::to_string(line_nr++).append(":").append("start"));
	   }

		in_file.close();
		parseOutputInfo();
		if (!CheckInput()) Input_error=true;
	} else {cout <<  "Inputfile " << name << " is not found. " << endl; Input_error=true; }

}
Input::~Input() {
}

bool Input::EvenSquareBrackets(const string& exp,vector<int> &open, vector<int> &close) const {
	return HasBalancedBrackets(exp, '[', ']', open, close);
}

bool Input::EvenBrackets(const string& exp,vector<int> &open, vector<int> &close) const {
	return HasBalancedBrackets(exp, '(', ')', open, close);
}

bool Input::ReadFile(const string& fname, string &In_buffer) const {
	ifstream this_file;
	bool success=true;
	const string resolved = ResolvePath(fname);
	this_file.open(resolved.c_str());
	std:: string In_line;
	if (this_file.is_open()) {
		while (this_file) {
			std::getline(this_file,In_line);
			NormalizeLine(In_line, false);
			if (In_line.empty()) continue;
			if (In_line.length()>2 && StartsWith(In_line, "//")) continue;
			In_buffer.append(In_line).append("#");
		}
		this_file.close();
		if (In_buffer.size()==0) {cout << "File " + resolved + " is empty " << endl; success=false; }
	} else {cout <<  "Inputfile " << resolved << " is not found. " << endl; success=false; }
	return success;
}


void Input::PrintList(const std::vector<std::string>& LIST) const {
	for (const std::string& item : LIST) {
		cout << item << " ; ";
	}
}

std::vector<std::string>& Input::split(const std::string& s, char delim, std::vector<std::string>&elems) const {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss,item,delim)) {
		item.erase(std::remove(item.begin(), item.end(), ' '), item.end());
		std::size_t pos = item.find("//");
		elems.push_back(item.substr(0,pos));
	}
	return elems;
}

bool Input:: TestNum(std::vector<std::string> &S, const string& c,int num_low, int num_high, int UptoStartNumber ) const {
	int n_starts=0;
	for (const std::string& entry : elems) {
		std::vector<std::string> set;
		split(entry, ':', set);
		if (set[1].substr(0,5)=="start") n_starts++;
		if (c==set[1] && n_starts<UptoStartNumber && !ContainsValue(S, set[2])) {
			S.push_back(set[2]);
		}
	}
	const int number = static_cast<int>(S.size());
	return number >= num_low && number <= num_high;
}

int Input:: GetNumStarts() const {
	int number=0;
	for (size_t i = 0; i < elems.size(); ++i) {
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1] == "start") number++;
		if (i == elems.size() - 1 && set[1].substr(0,5)!="start") {
			number++;
			//elems.push_back("start");
		}
	}
	return number;
}

bool Input:: InSet(const std::vector<std::string> &Standard, const string& keyword) const {
	return ContainsValue(Standard, keyword);
}
bool Input:: InSet(const std::vector<std::string> &Standard, int &pos, const string& keyword) const {
	return ContainsValue(Standard, keyword, &pos);
}

bool Input:: InSet(const vector<int> &Standard, int keyword) const {
	return ContainsValue(Standard, keyword);
}
bool Input:: InSet(const vector<int> &Standard, int &pos, int keyword) const {
	return ContainsValue(Standard, keyword, &pos);
}

// In->CheckParameters("keyword", name, start, KEYS, PARAMETERS)
bool Input::CheckParameters(const string& keyword, const string& name, int start, const std::vector<std::string>& Standard, ParameterStore& Input) const {
	bool success=true;
	bool prop_found;
	int length = elems.size();
	int S_length = Standard.size();
	int I_length;
	string parameter;
	int n_start=0;
	int n_found=0;
	int i=0;
	int j;
	std::vector<std::string> input_keys;
	std::vector<std::string> input_values;

	while (i<length && n_start<start){
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1]=="start") {
			n_start++;
			int k=0;
			int k_length=input_keys.size();
			while (k<k_length) { //remove doubles and keep last value; erase duplicates
				int l=k+1;
				int l_length=input_keys.size();
				while (l<l_length) {
					if (input_keys[k]==input_keys[l]) {
						input_values[k]=input_values[l];
						input_keys.erase(input_keys.begin()+l);
						if (n_start==1) cout <<"Warning: " << input_keys[k] << " found twice.... " << input_values[k] << " is used!" << endl;
						if (input_values.begin()+l != input_values.end()) input_values.erase(input_values.begin()+l);
						else input_values.erase(--input_values.end());
						l_length--;
						k_length--;
						l--;
					}
					l++;
				}
				k++;
			}
		}
		else {
			if (set[1] == keyword && set[2]== name) {
				parameter=set[3];
				j=0; prop_found=false;
				while (j<S_length && !prop_found) {
					if (Standard[j]==parameter) prop_found=true;
					j++;
				}
				if (!prop_found) {success=false; cout <<"In line " << set[0] << " "  << keyword << " property '" << parameter << "' is unknown. Select from: "<< endl;
					for (int k=0; k<S_length; k++) cout << Standard[k] << endl;
				} else {
					j=0; I_length = input_keys.size(); prop_found=false; n_found=0;
					while (j<I_length) {
						if (input_keys[j]==parameter) {prop_found = true; n_found++;}
						j++;
					}
					if (prop_found && n_found>1 && n_start==0) {success=false; cout <<n_start<<" "  << start << endl;  cout <<"In line " << set[0] << " " << keyword << " property '" << parameter << "' is already defined. "<< endl; }
					else {
						if (prop_found && n_found>1) {success=false; cout <<"After 'start' " << n_start << ", in line " << set[0] << " " << keyword << " property '" << parameter << "' is already defined. "<< endl; }
						else {input_keys.push_back(parameter); input_values.push_back(set[4]);
						}
					}
				}
			}
		}
		i++;
	}

	Input.clear();
	for (size_t idx = 0; idx < input_keys.size(); ++idx) {
		Input[input_keys[idx]] = input_values[idx];
	}

	return success;
}


bool Input:: LoadItems(const string& template_,std::vector<std::string> &Out_key, std::vector<std::string> &Out_name, std::vector<std::string> &Out_prop) const {
NAMICS_DBG("LoadItems in Input " << endl); Out_key.clear();
	Out_name.clear();
	Out_prop.clear();
	bool success=true;
	for (size_t i = 0; i < elems.size() ; i++) {
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1] == template_) {
			bool key_found=false;
			bool name_found=false;
			bool wild_monlist=false;
			bool wild_mollist=false;

			auto report_unknown_name = [&](const std::vector<std::string>& shown) {
				cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
				PrintList(shown);
			};
			auto validate_name = [&](const std::vector<std::string>& search, const std::vector<std::string>* shown = nullptr) {
				if (ContainsValue(search, set[3])) return true;
				report_unknown_name(shown == nullptr ? search : *shown);
				return false;
			};
			auto validate_singleton_name = [&](const std::vector<std::string>& singleton_list) {
				if (set[3]=="*" && !singleton_list.empty()) set[3]=singleton_list[0];
				return !singleton_list.empty() && set[3] == singleton_list[0] ? true : (report_unknown_name(singleton_list), false);
			};
			auto validate_name_or_wildcard = [&](const std::vector<std::string>& search, bool& wildcard) {
				if (set[3]=="*") {
					wildcard=true;
					return true;
				}
				return validate_name(search);
			};
			auto append_entries = [&](const std::vector<std::string>& names) {
				for (const std::string& name : names) {
					Out_key.push_back(set[2]);
					Out_name.push_back(name);
					Out_prop.push_back(set[4]);
				}
			};

			for (size_t j = 0 ; j < KEYS.size() ; j++) {
				if (KEYS[j] == set[2]) {
					key_found=true;
					switch (j-1) {
						case 0:
							name_found = validate_singleton_name(SysList);
							break;
						case 1:
							name_found = validate_name_or_wildcard(MolList, wild_mollist);
							break;
						case 2:
							name_found = validate_name_or_wildcard(MonList, wild_monlist);
							break;
						case 3:
							name_found = validate_singleton_name(LatList);
							break;
						case 4:
							name_found = validate_singleton_name(NewtonList);
							break;
						case 5:
							name_found = validate_name(OutputList);
							break;
						case 6:
							name_found=true;
							break;
						case 7:
							name_found = ContainsValue(StateList, set[3]);
							break;
						case 8:
							name_found = ContainsValue(ReactionList, set[3]);
							break;
						default:
							key_found=false;
						}

				}
			}
			if (!key_found) {cout<< "In line " << set[0] << " the keyword '" << set[2] << "' not recognized. Choose keywords from: " << endl;
				int length = KEYS.size();
				for (int k=1; k<length; k++) cout << KEYS[k] << endl;
				return false;
			}
			if (!name_found) {return false;}
			if (wild_mollist) {
				append_entries(MolList);
			}
			if (wild_monlist) {
				append_entries(MonList);
			}
			if (!(wild_monlist || wild_mollist)) {Out_key.push_back(set[2]); Out_name.push_back(set[3]); Out_prop.push_back(set[4]);}

		} //end temp
	} //end i;
	return success;
}

bool Input:: CheckInput(void) {
	bool success=true;
	int key_length=static_cast<int>(KEYS.size());
	bool last_start=false;

	if (elems.empty()) {cout << "inputfile is empty " << endl; success=false; }
	for (const std::string& entry : elems) {
		last_start=false;
		vector<std::string> set;
		split(entry,':',set);
		if (set.size() !=5) {if (set[1]!="start") {
			cout <<entry << endl;
			cout << " Line number " << set[0] << " does not contain 4 items" << endl; success=false; return false;}
		}
		if (set[1]=="start") last_start=true;
	}

	if (!last_start) {
		elems.push_back("0:start");
	}

	bool has_json_output = false;
	for (size_t i = 0; success && i < elems.size(); ++i) {
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1]=="output") {
			if (set[2] != "json") {
				cout << "Value for output extension '" << set[2] << "' not allowed. Only 'json' is supported." << endl;
				success=false;
			} else {
				has_json_output = true;
			}
		}
	}
	if (has_json_output && !InSet(KEYS, "json")) KEYS.push_back("json");
	key_length = static_cast<int>(KEYS.size());
	for (size_t i = 0; i < elems.size(); ++i) {
		vector<std::string> set;
		split(elems[i],':',set);
		const string& word=set[1];
		if (!InSet(KEYS, word)) {cout << word << " is not valid keyword in line " << set[0] << endl;
			cout << "select one of the following:" << endl;
			for( int k=0; k<key_length; k++) cout << KEYS[k] << endl;
			success=false;
		}
	}
	if (success) success=MakeLists(1);
	if (!OutputPathExists()) {
		cout << "Cannot access output folder '" << output_path << "'" << endl;
		success = false;
	}
	return success;
}

bool Input::MakeLists(int start) {
	bool success=true;
	SysList.clear();
	LatList.clear();
	NewtonList.clear();
	MonList.clear();
	MolList.clear();
	OutputList.clear();
	StateList.clear();
	ReactionList.clear();

	auto test_count = [&](std::vector<string>& list, const string& key, int low, int high, const string& error, bool fail_on_error = true) {
		const bool ok = TestNum(list, key, low, high, start);
		if (!ok && !error.empty()) cout << error << endl;
		if (!ok && fail_on_error) success=false;
		return ok;
	};

	test_count(SysList,"sys",0,1,"There can be no more than 1 'sys name' in the input");
	if (SysList.size()==0) SysList.push_back("NN");
	test_count(LatList,"lat",1,1,"There must be exactly one 'lat name' in the input");
	test_count(NewtonList,"newton",0,1,"There can be no more than 1 'newton name' in input");
	if (NewtonList.size()==0) NewtonList.push_back("NN");
	test_count(NewtonList,"newton",0,1,"There can be no more than 1 'newton name' in input");
	test_count(MonList,"mon",1,1000,"There must be at least one 'mon name' in input");
	test_count(StateList,"state",0,1000,"There can not be more than 1000 'state name's in input");
	test_count(ReactionList,"reaction",0,1000,"There can not be more than 1000 reaction name's in input");
	test_count(MolList,"mol",1,1000,"There must be at least one 'mol name' in input");
	test_count(OutputList,"output",1,1000,"No output defined! ", false);
	return success;
}

void Input::parseOutputInfo() {
	for (const string &line : elems) {
		vector<string> param;
		split(line, ':', param);
		if (param[1] != OUTPUT_INFO_KEY) {
			continue;
		}
		if (param[2] == "folder" && param[3] == "path") {
			output_path = NormalizeOutputPath(ResolvePath(param[4]));
		}
	}
}

const std::string& Input::GetOutputPath() const {
	return output_path;
}

std::string Input::ResolvePath(const std::string& path) const {
	if (path.empty()) return path;
	const std::filesystem::path candidate(path);
	if (candidate.is_absolute()) return candidate.lexically_normal().string();
	return (InputDirectoryFromName(name) / candidate).lexically_normal().string();
}

bool Input::OutputPathExists() const {
	return std::filesystem::is_directory(output_path);
}
