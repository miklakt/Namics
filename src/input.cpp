#include "input.h"

#include <filesystem>

namespace {

constexpr const char* OUTPUT_INFO_KEY = "out_info";
constexpr const char* DEFAULT_OUTPUT_PATH = "./output/";

const std::vector<std::string>& ProblemKeys() {
	static const std::vector<std::string> keys = {
		"sys", "mol", "mon", "lat", "newton", "output", OUTPUT_INFO_KEY, "state", "reaction", "sample", "json", "initial_guess"
	};
	return keys;
}

bool NeedsGuessOutputFile(const nlohmann::ordered_json& sys) {
	const auto write_guess = sys.find("write_initial_guess");
	if (write_guess == sys.end() || !write_guess->is_boolean() || !write_guess->get<bool>()) return false;
	const auto output_file = sys.find("guess_outputfile");
	return output_file == sys.end() || !output_file->is_string() || output_file->get<std::string>().empty();
}
} // namespace

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

Input::Input(const std::string& json_name)
	: json_path(json_name),
	  Input_error(false),
	  output_path(DEFAULT_OUTPUT_PATH),
	  starts(nlohmann::ordered_json::array()),
	  active_start(0) {
	std::ifstream input(json_path.c_str());
	if (!input.is_open()) {
		std::cout << "Inputfile " << json_path << " is not found. " << std::endl;
		Input_error = true;
		return;
	}

	nlohmann::ordered_json document;
	try {
		document = nlohmann::ordered_json::parse(input, nullptr, true, true);
	} catch (const std::exception& error) {
		std::cout << "Failed to parse JSON input file " << json_path << ": " << error.what() << std::endl;
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
	output_path = DEFAULT_OUTPUT_PATH;
	for (const auto& problem : starts) {
		if (!problem.is_object()) continue;
		const auto out_info = problem.find(OUTPUT_INFO_KEY);
		if (out_info == problem.end() || !out_info->is_object()) continue;
		const auto folder = out_info->find("folder");
		if (folder == out_info->end() || !folder->is_object()) continue;
		const auto path = folder->find("path");
		if (path == folder->end() || !path->is_string()) continue;
		output_path = ResolvePath(path->get<std::string>());
		if (output_path.empty()) output_path = DEFAULT_OUTPUT_PATH;
		else if (output_path.back() != '/') output_path.push_back('/');
	}
	if (!CheckInput()) Input_error = true;
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
	const auto& sections = (*this)[template_];
	if (!sections.is_object()) return out;

	auto fail = []() -> ParameterStore { return nlohmann::ordered_json(); };
	auto unknown_name = [](const std::string& name, const std::vector<std::string>& list) {
		std::cout << "Name '" << name << "' not recognised. Select from: " << std::endl;
		for (const std::string& item : list) std::cout << item << " ; ";
	};

	std::vector<std::string> names;
	std::vector<std::string> props;
	for (auto section_it = sections.begin(); section_it != sections.end(); ++section_it) {
		if (!section_it.value().is_object()) {
			std::cout << "json selector section '" << section_it.key() << "' must be an object." << std::endl;
			return fail();
		}
		for (auto entry_it = section_it.value().begin(); entry_it != section_it.value().end(); ++entry_it) {
			props.clear();
			if (entry_it.value().is_string()) props.push_back(entry_it.value().get<std::string>());
			else if (entry_it.value().is_array() && !entry_it.value().empty()) {
				for (const auto& prop : entry_it.value()) {
					if (!prop.is_string()) {
						std::cout << "json selector '" << section_it.key() << "." << entry_it.key()
						          << "' must only contain string properties." << std::endl;
						return fail();
					}
					props.push_back(prop.get<std::string>());
				}
			} else {
				std::cout << "json selector '" << section_it.key() << "." << entry_it.key()
				          << "' must be a string or a non-empty array of strings." << std::endl;
				return fail();
			}

			const std::vector<std::string>* list = nullptr;
			bool singleton = false;
			bool wildcard = false;
			if (section_it.key() == "sys") { list = &SysList; singleton = true; }
			else if (section_it.key() == "mol") { list = &MolList; wildcard = true; }
			else if (section_it.key() == "mon") { list = &MonList; wildcard = true; }
			else if (section_it.key() == "lat") { list = &LatList; singleton = true; }
			else if (section_it.key() == "newton") { list = &NewtonList; singleton = true; }
			else if (section_it.key() == "output") list = &OutputList;
			else if (section_it.key() == OUTPUT_INFO_KEY) names = {entry_it.key()};
			else if (section_it.key() == "state") list = &StateList;
			else if (section_it.key() == "reaction") list = &ReactionList;
			else {
				std::cout << "The selector section '" << section_it.key() << "' is not recognized." << std::endl;
				return fail();
			}

			if (section_it.key() != OUTPUT_INFO_KEY) {
				names.clear();
				std::string name_ = entry_it.key();
				if (singleton && name_ == "*" && !list->empty()) name_ = (*list)[0];
				if (wildcard && name_ == "*") names = *list;
				else if (ContainsValue(*list, name_)) names.push_back(name_);
				else {
					unknown_name(entry_it.key(), *list);
					return fail();
				}
			}

			for (const std::string& name_ : names) {
				for (const std::string& prop : props) {
					if (prop.find('(') != std::string::npos || prop.find(')') != std::string::npos) {
						std::cout << "Indexed json output selectors like '" << prop << "' are no longer supported." << std::endl;
						return fail();
					}
					out.push_back({{"key", section_it.key()}, {"name", name_}, {"prop", prop}});
				}
			}
		}
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
				if (!it.value().is_object()) {
					std::cout << "'json' in problem " << i + 1 << " must be an object." << std::endl;
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
			if (it.key() == "sys") {
				for (auto sys_it = it.value().begin(); sys_it != it.value().end(); ++sys_it) {
					if (!sys_it.value().is_object()) {
						std::cout << "sys." << sys_it.key() << " in problem " << i + 1 << " must be an object." << std::endl;
						success = false;
						continue;
					}
					if (!NeedsGuessOutputFile(sys_it.value())) continue;
					std::cout << "When 'write_initial_guess' is true in problem " << i + 1
					          << ", provide 'guess_outputfile'." << std::endl;
					success = false;
				}
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
	if (!std::filesystem::is_directory(output_path)) {
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

std::string Input::ResolvePath(const std::string& path) const {
	if (path.empty()) return path;
	const std::filesystem::path candidate(path);
	if (candidate.is_absolute()) return candidate.lexically_normal().string();
	const std::filesystem::path base(json_path);
	const std::filesystem::path parent = base.has_parent_path() ? base.parent_path() : std::filesystem::path(".");
	return (parent / candidate).lexically_normal().string();
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
