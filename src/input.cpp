#include "input.h"

#include <filesystem>

namespace {

constexpr const char* OUTPUT_INFO_KEY = "out_info";
constexpr const char* DEFAULT_OUTPUT_PATH = "./output/";

const std::vector<std::string>& ProblemKeys() {
	static const std::vector<std::string> keys = {
		"sys", "mol", "mon", "lat", "newton", "output", OUTPUT_INFO_KEY, "state", "reaction", "json", "initial_guess"
	};
	return keys;
}

std::string NormalizeOutputPath(std::string path) {
	if (path.empty()) return DEFAULT_OUTPUT_PATH;
	if (path.back() != '/') path.push_back('/');
	return path;
}

bool HasBalancedBrackets(
	const std::string& expression,
	char open_bracket,
	char close_bracket,
	std::vector<int>& open_positions,
	std::vector<int>& close_positions) {
	std::vector<char> stack;
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

Input::Input(const std::string& json_name)
	: json_path(json_name),
	  Input_error(false),
	  output_path(DEFAULT_OUTPUT_PATH),
	  starts(nlohmann::ordered_json::array()),
	  active_start(0) {
	const std::string load_name = json_path;
	std::ifstream input(load_name.c_str());
	if (!input.is_open()) {
		std::cout << "Inputfile " << load_name << " is not found. " << std::endl;
		Input_error = true;
		return;
	}

	nlohmann::ordered_json document;
	try {
		document = nlohmann::ordered_json::parse(input, nullptr, true, true);
	} catch (const std::exception& error) {
		std::cout << "Failed to parse JSON input file " << load_name << ": " << error.what() << std::endl;
		Input_error = true;
		return;
	}

	if (document.is_array()) {
		starts = std::move(document);
	} else if (document.is_object()) {
		const auto problems = document.find("problems");
		if (problems == document.end()) {
			starts = nlohmann::ordered_json::array({std::move(document)});
		} else if (problems->is_array()) {
			starts = *problems;
		} else {
			starts = nlohmann::ordered_json::array({*problems});
		}
	} else {
		std::cout << "JSON input root must be a problem object, a problem array, or an object containing 'problems'." << std::endl;
		Input_error = true;
		return;
	}
	UpdateOutputPath();
	if (!CheckInput()) Input_error = true;
}

Input::~Input() {
}

bool Input::EvenSquareBrackets(const std::string& exp, std::vector<int>& open, std::vector<int>& close) const {
	return HasBalancedBrackets(exp, '[', ']', open, close);
}

bool Input::EvenBrackets(const std::string& exp, std::vector<int>& open, std::vector<int>& close) const {
	return HasBalancedBrackets(exp, '(', ')', open, close);
}

std::vector<std::string>& Input::split(const std::string& s, char delim, std::vector<std::string>& elems) const {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		item.erase(std::remove(item.begin(), item.end(), ' '), item.end());
		const std::size_t pos = item.find("//");
		elems.push_back(item.substr(0, pos));
	}
	return elems;
}

int Input::GetNumStarts() const {
	return starts.is_array() ? static_cast<int>(starts.size()) : 0;
}

const ParameterStore& Input::Parameters(const std::string& keyword, const std::string& name_, int start) const {
	static const ParameterStore empty = nlohmann::ordered_json::object();
	const auto& start_object = Start(start);
	if (!start_object.is_object()) return empty;
	const auto keyword_it = start_object.find(keyword);
	if (keyword_it == start_object.end() || !keyword_it->is_object()) return empty;
	const auto name_it = keyword_it->find(name_);
	if (name_it == keyword_it->end() || !name_it->is_object()) return empty;
	return name_it.value();
}

ParameterStore Input::LoadItems(const std::string& template_) const {
NAMICS_DBG("LoadItems in Input " << std::endl);
	ParameterStore out = nlohmann::ordered_json::array();
	const auto& items = (*this)[template_];
	if (!items.is_array()) return out;
	for (const auto& item : items) {
		std::string key = item.value("key", "");
		std::string name_ = item.value("name", "");
		const std::string prop = item.value("prop", "");
		if (prop.find('(') != std::string::npos || prop.find(')') != std::string::npos) {
			std::cout << "Indexed json output selectors like '" << prop << "' are no longer supported." << std::endl;
			return nlohmann::ordered_json();
		}
		bool key_found = false;
		bool name_found = false;
		bool wild_monlist = false;
		bool wild_mollist = false;

		auto report_unknown_name = [&](const std::vector<std::string>& shown) {
			std::cout << "Name '" << name_ << "' not recognised. Select from: " << std::endl;
			for (const std::string& item : shown) std::cout << item << " ; ";
		};
		auto validate_name = [&](const std::vector<std::string>& search, const std::vector<std::string>* shown = nullptr) {
			if (ContainsValue(search, name_)) return true;
			report_unknown_name(shown == nullptr ? search : *shown);
			return false;
		};
		auto validate_singleton_name = [&](const std::vector<std::string>& singleton_list) {
			if (name_ == "*" && !singleton_list.empty()) name_ = singleton_list[0];
			if (!singleton_list.empty() && name_ == singleton_list[0]) return true;
			report_unknown_name(singleton_list);
			return false;
		};
		auto validate_name_or_wildcard = [&](const std::vector<std::string>& search, bool& wildcard) {
			if (name_ == "*") {
				wildcard = true;
				return true;
			}
			return validate_name(search);
		};
		auto append_entry = [&](const std::string& expanded_name) {
			out.push_back({
				{"key", key},
				{"name", expanded_name},
				{"prop", prop}
			});
		};

		if (key == "sys") {
			key_found = true;
			name_found = validate_singleton_name(SysList);
		} else if (key == "mol") {
			key_found = true;
			name_found = validate_name_or_wildcard(MolList, wild_mollist);
		} else if (key == "mon") {
			key_found = true;
			name_found = validate_name_or_wildcard(MonList, wild_monlist);
		} else if (key == "lat") {
			key_found = true;
			name_found = validate_singleton_name(LatList);
		} else if (key == "newton") {
			key_found = true;
			name_found = validate_singleton_name(NewtonList);
		} else if (key == "output") {
			key_found = true;
			name_found = validate_name(OutputList);
		} else if (key == OUTPUT_INFO_KEY) {
			key_found = true;
			name_found = true;
		} else if (key == "state") {
			key_found = true;
			name_found = ContainsValue(StateList, name_);
		} else if (key == "reaction") {
			key_found = true;
			name_found = ContainsValue(ReactionList, name_);
		}

		if (!key_found) {
			std::cout << "The keyword '" << key << "' not recognized. Choose keywords from: " << std::endl;
			for (const std::string& item_name : ProblemKeys()) std::cout << item_name << std::endl;
			return nlohmann::ordered_json();
		}
		if (!name_found) return nlohmann::ordered_json();
		if (wild_mollist) {
			for (const std::string& mol_name : MolList) append_entry(mol_name);
		}
		if (wild_monlist) {
			for (const std::string& mon_name : MonList) append_entry(mon_name);
		}
		if (!(wild_monlist || wild_mollist)) append_entry(name_);
	}
	return out;
}

bool Input::CheckInput() {
	bool success = true;
	if (!starts.is_array() || starts.empty()) {
		std::cout << "inputfile is empty " << std::endl;
		success = false;
	}

	for (size_t i = 0; i < starts.size(); ++i) {
		const auto& problem = starts[i];
		if (!problem.is_object()) {
			std::cout << "Problem " << i + 1 << " is not a JSON object." << std::endl;
			success = false;
			continue;
		}
		for (auto it = problem.begin(); it != problem.end(); ++it) {
			if (!ContainsValue(ProblemKeys(), it.key())) {
				std::cout << it.key() << " is not a valid problem key in problem " << i + 1 << std::endl;
				std::cout << "select one of the following:" << std::endl;
				for (const std::string& item : ProblemKeys()) std::cout << item << std::endl;
				success = false;
				continue;
			}
			if (it.key() == "json") {
				if (!it.value().is_array()) {
					std::cout << "'json' in problem " << i + 1 << " must be an array." << std::endl;
					success = false;
				}
				continue;
			}
			if (it.key() == "initial_guess") {
				if (!it.value().is_object()) {
					std::cout << "'initial_guess' in problem " << i + 1 << " must be an object." << std::endl;
					success = false;
				}
				continue;
			}
			if (!it.value().is_object()) {
				std::cout << "'" << it.key() << "' in problem " << i + 1 << " must be an object." << std::endl;
				success = false;
				continue;
			}
			if (it.key() != "output") continue;
			for (auto output_it = it.value().begin(); output_it != it.value().end(); ++output_it) {
				if (output_it.key() != "json") {
					std::cout << "Value for output extension '" << output_it.key() << "' not allowed. Only 'json' is supported." << std::endl;
					success = false;
				}
				if (!output_it.value().is_object()) {
					std::cout << "output." << output_it.key() << " in problem " << i + 1 << " must be an object." << std::endl;
					success = false;
				}
			}
		}
	}

	if (success) success = MakeLists(1);
	if (!OutputPathExists()) {
		std::cout << "Cannot access output folder '" << output_path << "'" << std::endl;
		success = false;
	}
	return success;
}

bool Input::MakeLists(int start) {
	bool success = true;
	active_start = std::max(0, start - 1);
	SysList.clear();
	LatList.clear();
	NewtonList.clear();
	MonList.clear();
	MolList.clear();
	OutputList.clear();
	StateList.clear();
	ReactionList.clear();

	const auto collect_names = [&](const std::string& key, std::vector<std::string>& list) {
		const auto& section = (*this)[key];
		if (!section.is_object()) return;
		for (auto it = section.begin(); it != section.end(); ++it) list.push_back(it.key());
	};

	collect_names("sys", SysList);
	collect_names("lat", LatList);
	collect_names("newton", NewtonList);
	collect_names("mon", MonList);
	collect_names("mol", MolList);
	collect_names("output", OutputList);
	collect_names("state", StateList);
	collect_names("reaction", ReactionList);

	auto test_count = [&](std::vector<std::string>& list, int low, int high, const std::string& error, bool fail_on_error = true) {
		const bool ok = static_cast<int>(list.size()) >= low && static_cast<int>(list.size()) <= high;
		if (!ok && !error.empty()) std::cout << error << std::endl;
		if (!ok && fail_on_error) success = false;
		return ok;
	};

	test_count(SysList, 0, 1, "There can be no more than 1 'sys name' in the input");
	if (SysList.empty()) SysList.push_back("NN");
	test_count(LatList, 1, 1, "There must be exactly one 'lat name' in the input");
	test_count(NewtonList, 0, 1, "There can be no more than 1 'newton name' in input");
	if (NewtonList.empty()) NewtonList.push_back("NN");
	test_count(MonList, 1, 1000, "There must be at least one 'mon name' in input");
	test_count(StateList, 0, 1000, "There can not be more than 1000 'state name's in input");
	test_count(ReactionList, 0, 1000, "There can not be more than 1000 reaction name's in input");
	test_count(MolList, 1, 1000, "There must be at least one 'mol name' in input");
	test_count(OutputList, 1, 1000, "No output defined! ", false);
	return success;
}

void Input::UpdateOutputPath() {
	output_path = DEFAULT_OUTPUT_PATH;
	for (const auto& problem : starts) {
		if (!problem.is_object()) continue;
		const auto out_info = problem.find(OUTPUT_INFO_KEY);
		if (out_info == problem.end() || !out_info->is_object()) continue;
		const auto folder = out_info->find("folder");
		if (folder == out_info->end() || !folder->is_object()) continue;
		const auto path = folder->find("path");
		if (path == folder->end() || !path->is_string()) continue;
		output_path = NormalizeOutputPath(ResolvePath(path->get<std::string>()));
	}
}

const std::string& Input::GetOutputPath() const {
	return output_path;
}

std::string Input::ResolvePath(const std::string& path) const {
	if (path.empty()) return path;
	const std::filesystem::path candidate(path);
	if (candidate.is_absolute()) return candidate.lexically_normal().string();
	const std::filesystem::path base(json_path);
	const std::filesystem::path parent = base.has_parent_path() ? base.parent_path() : std::filesystem::path(".");
	return (parent / candidate).lexically_normal().string();
}

bool Input::OutputPathExists() const {
	return std::filesystem::is_directory(output_path);
}

const ParameterStore& Input::operator[](const std::string& key) const {
	const auto& start_object = Start(active_start + 1);
	static const ParameterStore empty = nlohmann::ordered_json::object();
	if (!start_object.is_object()) return empty;
	const auto it = start_object.find(key);
	return it == start_object.end() ? empty : it.value();
}

const ParameterStore& Input::Start(int start) const {
	static const ParameterStore empty = nlohmann::ordered_json::object();
	if (!starts.is_array() || start < 1 || start > static_cast<int>(starts.size())) return empty;
	return starts[static_cast<size_t>(start - 1)];
}
