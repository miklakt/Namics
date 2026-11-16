#include "input.h"

#include <filesystem>

namespace {

const std::vector<std::string>& ProblemKeys() {
	static const std::vector<std::string> keys = {
		"sys", "mol", "mon", "lat", "newton", "state", "reaction", "sample", "initial_guess"
	};
	return keys;
}

bool NeedsGuessOutputFile(const nlohmann::ordered_json& sys) {
	const auto write_guess = sys.find("write_initial_guess");
	if (write_guess == sys.end() || !write_guess->is_boolean() || !write_guess->get<bool>()) return false;
	const auto output_file = sys.find("guess_outputfile");
	return output_file == sys.end() || !output_file->is_string() || output_file->get<std::string>().empty();
}

std::string OutputFile(const std::string& input_path, const std::string& configured) {
	const std::string suffix = ".output.json";
	std::string file = configured;
	if (file.empty()) {
		file = input_path;
		for (const std::string ext : {".input.json", ".json", ".in"}) {
			if (file.size() > ext.size() && file.compare(file.size() - ext.size(), ext.size(), ext) == 0) {
				file.resize(file.size() - ext.size());
				break;
			}
		}
	}
	if (file.size() < suffix.size() || file.compare(file.size() - suffix.size(), suffix.size(), suffix) != 0) file += suffix;
	return file;
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

Input::Input(const std::string& path_, ParameterStore document)
	: path(path_),
	  Input_error(false),
	  output(nlohmann::ordered_json::object()),
	  starts(std::move(document["problems"])),
	  active_start(0) {
	for (auto it = document.begin(); it != document.end(); ++it) {
		if (it.key() == "problems") continue;
		if (it.key() == "output" && it.value().is_object()) {
			output = std::move(it.value());
			continue;
		}
		std::cout << "Input root key '" << it.key() << "' is not allowed; expected 'output' (object) and 'problems'." << std::endl;
		Input_error = true;
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

ParameterStore Input::Selectors() const {
	ParameterStore out = nlohmann::ordered_json::array();
	const auto sections = output.find("select");
	if (sections == output.end()) return out;
	if (!sections->is_object()) {
		std::cout << "output.select must be an object." << std::endl;
		return nlohmann::ordered_json();
	}

	auto fail = []() -> ParameterStore { return nlohmann::ordered_json(); };
	std::vector<std::string> names;
	std::vector<std::string> props;
	for (auto section_it = sections->begin(); section_it != sections->end(); ++section_it) {
		if (!section_it.value().is_object()) {
			std::cout << "output.select." << section_it.key() << " must be an object." << std::endl;
			return fail();
		}
		const std::vector<std::string>* list = nullptr;
		bool singleton = false;
		if (section_it.key() == "sys") { list = &SysList; singleton = true; }
		else if (section_it.key() == "mol") list = &MolList;
		else if (section_it.key() == "mon") list = &MonList;
		else if (section_it.key() == "lat") { list = &LatList; singleton = true; }
		else if (section_it.key() == "newton") { list = &NewtonList; singleton = true; }
		else if (section_it.key() == "state") list = &StateList;
		else if (section_it.key() == "reaction") list = &ReactionList;
		else {
			std::cout << "The selector section '" << section_it.key() << "' is not recognized." << std::endl;
			return fail();
		}
		for (auto entry_it = section_it.value().begin(); entry_it != section_it.value().end(); ++entry_it) {
			const std::string label = "output.select." + section_it.key() + "." + entry_it.key();
			props.clear();
			if (entry_it.value().is_string()) props.push_back(entry_it.value().get<std::string>());
			else if (entry_it.value().is_array() && !entry_it.value().empty()) {
				for (const auto& prop : entry_it.value()) {
					if (!prop.is_string()) {
						std::cout << label << " must only contain string properties." << std::endl;
						return fail();
					}
					props.push_back(prop.get<std::string>());
				}
			} else {
				std::cout << label << " must be a string or a non-empty array of strings." << std::endl;
				return fail();
			}

			names.clear();
			const std::string& name_ = entry_it.key();
			if (singleton) names.push_back((*list)[0]);
			else if (name_ == "*") names = *list;
			else if (ContainsValue(*list, name_)) names.push_back(name_);
			else {
				std::cout << "Name '" << name_ << "' in " << label << " not recognised. Select from: " << std::endl;
				for (const std::string& item : *list) std::cout << item << " ; ";
				std::cout << std::endl;
				return fail();
			}

			for (const std::string& name : names) {
				for (const std::string& prop : props) {
					if (prop.find('(') != std::string::npos || prop.find(')') != std::string::npos) {
						std::cout << "Indexed json output selectors like '" << prop << "' are no longer supported." << std::endl;
						return fail();
					}
					out.push_back({{"key", section_it.key()}, {"name", name}, {"prop", prop}});
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
			if (!it.value().is_object()) {
				std::cout << "'" << it.key() << "' in problem " << i + 1 << " must be an object." << std::endl;
				success = false;
				continue;
			}
			if (it.key() != "sys") continue;
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
		}
	}

	static const std::vector<std::string> output_keys = {"file", "append", "write_bounds", "select"};
	for (auto it = output.begin(); it != output.end(); ++it) {
		if (ContainsValue(output_keys, it.key())) continue;
		std::cout << "output property '" << it.key() << "' is unknown. Select from: file, append, write_bounds, select" << std::endl;
		success = false;
	}
	try {
		output_file = OutputFile(path, ResolvePath(output.value("file", std::string{})));
		(void)output.value("append", false);
		(void)output.value("write_bounds", false);
	} catch (const nlohmann::json::exception& error) {
		std::cout << "Invalid json type in output: " << error.what() << std::endl;
		success = false;
	}

	return success && MakeLists(1);
}

bool Input::MakeLists(int start) {
	bool success = true;
	active_start = std::max(0, start - 1);
	SysList.clear();
	LatList.clear();
	NewtonList.clear();
	MonList.clear();
	MolList.clear();
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
	collect_names("state", StateList);
	collect_names("reaction", ReactionList);

	auto test_count = [&](std::vector<std::string>& list, int low, int high, const std::string& error) {
		if (static_cast<int>(list.size()) >= low && static_cast<int>(list.size()) <= high) return;
		std::cout << error << std::endl;
		success = false;
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
	return success;
}

std::string Input::ResolvePath(const std::string& relative) const {
	if (relative.empty()) return relative;
	const std::filesystem::path candidate(relative);
	if (candidate.is_absolute()) return candidate.lexically_normal().string();
	const std::filesystem::path base(path);
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
