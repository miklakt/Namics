#include "input_preprocessor.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>

namespace {

using json = nlohmann::ordered_json;
using Coordinates = std::vector<std::vector<int>>;

struct ProblemShape {
	int gradients = 1;
	std::array<int, 3> layers = {0, 1, 1};
	std::array<std::string, 6> bc = {"mirror", "mirror", "mirror", "mirror", "mirror", "mirror"};
};

constexpr std::array<const char*, 3> MASK_FILE_KEYS = {"file", "filename", "path"};
constexpr const char* MASK_COORDINATES_KEY = "coordinates";
constexpr std::array<const char*, 17> LEGACY_RANGE_KEYWORDS = {
	"firstlayer", "lastlayer", "var_pos",
	"firstlayer_x", "firstlayer_y", "firstlayer_z",
	"lastlayer_x", "lastlayer_y", "lastlayer_z",
	"lowerbound", "upperbound",
	"lowerbound_x", "lowerbound_y", "lowerbound_z",
	"upperbound_x", "upperbound_y", "upperbound_z"
};

bool StartsWith(const std::string& value, const std::string& prefix) {
	return value.rfind(prefix, 0) == 0;
}

bool EndsWith(const std::string& value, const std::string& suffix) {
	return value.size() >= suffix.size() &&
		value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void NormalizeLegacyLine(std::string& line) {
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
	line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
}

std::string ToLowerCopy(std::string value) {
	std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
		return static_cast<char>(std::tolower(c));
	});
	return value;
}

template <typename T>
bool ParseStrict(const std::string& s, T& value) {
	if (s.empty()) return false;
	std::istringstream stream(s);
	stream >> std::noskipws >> value;
	return !stream.fail() && stream.eof();
}

template <>
bool ParseStrict<bool>(const std::string& s, bool& value) {
	const std::string lowered = ToLowerCopy(s);
	if (lowered == "true") {
		value = true;
		return true;
	}
	if (lowered == "false") {
		value = false;
		return true;
	}
	return false;
}

bool IsEmptyOrComment(const std::string& line) {
	return line.empty() || StartsWith(line, "//");
}

std::string TrimCopy(std::string value) {
	const auto not_space = [](unsigned char c) { return !std::isspace(c); };
	value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
	value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
	return value;
}

std::vector<std::string> SplitLegacy(const std::string& line, char delim) {
	std::vector<std::string> parts;
	std::stringstream ss(line);
	std::string item;
	while (std::getline(ss, item, delim)) {
		item.erase(std::remove(item.begin(), item.end(), ' '), item.end());
		const std::size_t pos = item.find("//");
		parts.push_back(item.substr(0, pos));
	}
	return parts;
}

json ParseLegacyValue(const std::string& value) {
	bool bool_value = false;
	long long int_value = 0;
	double real_value = 0;
	if (ParseStrict(value, bool_value)) return bool_value;
	if (ParseStrict(value, int_value)) return int_value;
	if (ParseStrict(value, real_value)) return real_value;
	return value;
}

std::filesystem::path ResolveRelativeTo(const std::string& base_path, const std::string& other_path) {
	const std::filesystem::path candidate(other_path);
	if (candidate.is_absolute()) return candidate.lexically_normal();
	const std::filesystem::path base(base_path);
	const std::filesystem::path parent = base.has_parent_path() ? base.parent_path() : std::filesystem::path(".");
	return (parent / candidate).lexically_normal();
}

bool ReadJsonFile(const std::string& path, json& document) {
	std::ifstream input(path.c_str());
	if (!input.is_open()) {
		std::cout << "Inputfile " << path << " is not found. " << std::endl;
		return false;
	}
	try {
		document = json::parse(input, nullptr, true, true);
	} catch (const std::exception& error) {
		std::cout << "Failed to parse JSON input file " << path << ": " << error.what() << std::endl;
		return false;
	}
	return true;
}

bool LooksLikeJsonInput(const std::string& path) {
	std::ifstream in(path.c_str());
	if (!in.is_open()) return false;
	char ch = '\0';
	while (in.get(ch)) {
		if (std::isspace(static_cast<unsigned char>(ch))) continue;
		if (ch == '/' && in.peek() == '/') {
			in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		return ch == '{' || ch == '[';
	}
	return false;
}

bool ReadLegacyEntries(const std::string& path, std::vector<std::string>& entries) {
	std::ifstream in(path.c_str());
	if (!in.is_open()) {
		std::cout << "Inputfile " << path << " is not found. " << std::endl;
		return false;
	}
	std::string line;
	while (std::getline(in, line)) {
		NormalizeLegacyLine(line);
		if (IsEmptyOrComment(line)) continue;
		if (StartsWith(line, "include:")) {
			std::string include_spec = line.substr(8);
			if (EndsWith(include_spec, "::")) include_spec.resize(include_spec.size() - 2);
			const std::string include_path = ResolveRelativeTo(path, include_spec).string();
			std::ifstream include_file(include_path.c_str());
			if (!include_file.is_open()) {
				std::cout << "'include : " << include_path << "' not working because file is not found " << std::endl;
				return false;
			}
			std::string include_line;
			while (std::getline(include_file, include_line)) {
				NormalizeLegacyLine(include_line);
				if (!IsEmptyOrComment(include_line)) entries.push_back(include_line);
			}
			continue;
		}
		entries.push_back(line);
	}
	if (!entries.empty() && entries.back() != "start") entries.push_back("start");
	return true;
}

void AppendJsonSelector(json& problem, const std::string& section, const std::string& entry, const std::string& property) {
	auto& selectors = problem["json"];
	if (!selectors.is_object()) selectors = json::object();
	auto& section_selectors = selectors[section];
	if (!section_selectors.is_object()) section_selectors = json::object();
	auto& properties = section_selectors[entry];
	if (properties.is_null()) {
		properties = property;
		return;
	}
	if (properties.is_string()) {
		properties = json::array({properties.get<std::string>(), property});
		return;
	}
	if (!properties.is_array()) properties = json::array();
	properties.push_back(property);
}

bool ConvertLegacyEntries(const std::vector<std::string>& entries, json& problems) {
	json current = json::object();
	for (const std::string& entry : entries) {
		if (entry == "start") {
			problems.push_back(current);
			continue;
		}
		const std::vector<std::string> parts = SplitLegacy(entry, ':');
		if (parts.empty()) continue;
		if (parts[0] == "json") {
			if (parts.size() != 4) return false;
			AppendJsonSelector(current, parts[1], parts[2], parts[3]);
			continue;
		}
		if (parts.size() != 4) return false;
		auto& section = current[parts[0]];
		if (!section.is_object()) section = json::object();
		auto& object = section[parts[1]];
		if (!object.is_object()) object = json::object();
		object[parts[2]] = ParseLegacyValue(parts[3]);
	}
	return true;
}

std::filesystem::path EmittedJsonPath(const std::string& source_path) {
	const std::filesystem::path source(source_path);
	if (EndsWith(source.filename().string(), ".input.json")) return source;
	const std::filesystem::path parent = source.has_parent_path() ? source.parent_path() : std::filesystem::path(".");
	if (source.has_stem()) return parent / (source.stem().string() + ".input.json");
	std::filesystem::path emitted = source.filename();
	emitted += ".input.json";
	return parent / emitted;
}

bool ReadSanitizedTextFile(const std::string& path, std::vector<std::string>& tokens) {
	std::ifstream input(path.c_str());
	if (!input.is_open()) {
		std::cout << "Inputfile " << path << " is not found. " << std::endl;
		return false;
	}
	tokens.clear();
	std::string line;
	while (std::getline(input, line)) {
		NormalizeLegacyLine(line);
		if (IsEmptyOrComment(line)) continue;
		std::stringstream ss(line);
		std::string token;
		while (std::getline(ss, token, '#')) {
			if (!token.empty()) tokens.push_back(token);
		}
	}
	return true;
}

template <size_t N>
const json* FindMember(const json& object, const std::array<const char*, N>& keys) {
	if (!object.is_object()) return nullptr;
	for (const char* key : keys) {
		const auto it = object.find(key);
		if (it != object.end()) return &it.value();
	}
	return nullptr;
}

bool ParseCoordinateList(const json& value, int dimensions, Coordinates& coordinates) {
	if (!value.is_array()) return false;
	try {
		coordinates = value.get<Coordinates>();
	} catch (const json::exception&) {
		coordinates.clear();
		return false;
	}
	for (const auto& point : coordinates) {
		if (static_cast<int>(point.size()) == dimensions) continue;
		coordinates.clear();
		return false;
	}
	return true;
}

bool LoadProblemShape(const json& problem, ProblemShape& shape) {
	const auto lat_section = problem.find("lat");
	if (lat_section == problem.end() || !lat_section->is_object() || lat_section->empty()) return false;
	const auto parameters_it = lat_section->begin();
	if (!parameters_it.value().is_object()) return false;
	const json& parameters = parameters_it.value();
	try {
		shape.gradients = parameters.value("gradients", 1);
		if (shape.gradients < 1 || shape.gradients > 3) return false;
		shape.layers[0] = shape.gradients == 1 ? parameters.value("n_layers", -1) : parameters.value("n_layers_x", -1);
		shape.layers[1] = shape.gradients >= 2 ? parameters.value("n_layers_y", -1) : 1;
		shape.layers[2] = shape.gradients >= 3 ? parameters.value("n_layers_z", -1) : 1;
		if (shape.layers[0] < 1 || (shape.gradients >= 2 && shape.layers[1] < 1) || (shape.gradients >= 3 && shape.layers[2] < 1)) {
			return false;
		}
		if (shape.gradients == 1) {
			shape.bc[0] = parameters.value("lowerbound", std::string{"mirror"});
			shape.bc[3] = parameters.value("upperbound", std::string{"mirror"});
		} else {
			shape.bc[0] = parameters.value("lowerbound_x", std::string{"mirror"});
			shape.bc[3] = parameters.value("upperbound_x", std::string{"mirror"});
			shape.bc[1] = parameters.value("lowerbound_y", std::string{"mirror"});
			shape.bc[4] = parameters.value("upperbound_y", std::string{"mirror"});
			if (shape.gradients == 3) {
				shape.bc[2] = parameters.value("lowerbound_z", std::string{"mirror"});
				shape.bc[5] = parameters.value("upperbound_z", std::string{"mirror"});
			}
		}
	} catch (const json::exception&) {
		return false;
	}
	return true;
}

bool ParseLegacyMaskFile(const std::string& path, const ProblemShape& shape, Coordinates& coordinates) {
	std::vector<std::string> tokens;
	if (!ReadSanitizedTextFile(path, tokens)) return false;
	coordinates.clear();
	for (const std::string& token : tokens) {
		const std::vector<std::string> parts = SplitLegacy(token, ',');
		if (static_cast<int>(parts.size()) != shape.gradients) return false;
		std::vector<int> point;
		point.reserve(static_cast<size_t>(shape.gradients));
		for (const std::string& part : parts) {
			int parsed = 0;
			if (!ParseStrict(part, parsed)) return false;
			point.push_back(parsed);
		}
		coordinates.push_back(std::move(point));
	}
	return true;
}

bool IsLegacyRangeText(const std::string& value) {
	if (value.find(';') != std::string::npos || value.find(',') != std::string::npos) return true;
	const std::string lowered = ToLowerCopy(TrimCopy(value));
	return std::find(LEGACY_RANGE_KEYWORDS.begin(), LEGACY_RANGE_KEYWORDS.end(), lowered) != LEGACY_RANGE_KEYWORDS.end();
}

bool ExpandLegacyRangeAlias(const std::string& alias,
                            const std::string& kind,
                            const ProblemShape& shape,
                            std::vector<int>& first,
                            std::vector<int>& last) {
	first.assign(static_cast<size_t>(shape.gradients), 1);
	last.assign(static_cast<size_t>(shape.gradients), 1);
	for (int axis = 0; axis < shape.gradients; ++axis) last[static_cast<size_t>(axis)] = shape.layers[axis];
	const std::string lowered = ToLowerCopy(alias);
	if (kind == "pinned") {
		if (lowered == "firstlayer" || lowered == "firstlayer_x") last[0] = first[0] = 1;
		else if (lowered == "lastlayer" || lowered == "lastlayer_x") first[0] = last[0] = shape.layers[0];
		else if ((lowered == "firstlayer_y" || lowered == "lastlayer_y") && shape.gradients >= 2) {
			first[1] = last[1] = lowered == "firstlayer_y" ? 1 : shape.layers[1];
		} else if ((lowered == "firstlayer_z" || lowered == "lastlayer_z") && shape.gradients >= 3) {
			first[2] = last[2] = lowered == "firstlayer_z" ? 1 : shape.layers[2];
		} else {
			return false;
		}
		return true;
	}
	auto set_surface = [&](int axis, bool upper) {
		first[axis] = last[axis] = upper ? shape.layers[axis] + 1 : 0;
		return shape.bc[static_cast<size_t>(axis + (upper ? 3 : 0))] == "surface";
	};
	if (lowered == "lowerbound" || lowered == "lowerbound_x") return set_surface(0, false);
	if (lowered == "upperbound" || lowered == "upperbound_x") return set_surface(0, true);
	if ((lowered == "lowerbound_y" || lowered == "upperbound_y") && shape.gradients >= 2) return set_surface(1, lowered == "upperbound_y");
	if ((lowered == "lowerbound_z" || lowered == "upperbound_z") && shape.gradients >= 3) return set_surface(2, lowered == "upperbound_z");
	return false;
}

bool ParseLegacyRange(const std::string& value,
                      const std::string& kind,
                      const ProblemShape& shape,
                      int var_pos,
                      Coordinates& coordinates) {
	std::string compact = TrimCopy(value);
	compact.erase(std::remove(compact.begin(), compact.end(), ' '), compact.end());
	while (!compact.empty() && compact.back() == ';') compact.pop_back();
	std::vector<int> first;
	std::vector<int> last;
	if (!ExpandLegacyRangeAlias(compact, kind, shape, first, last)) {
		const std::vector<std::string> endpoints = SplitLegacy(compact, ';');
		if (endpoints.size() != 2) return false;
		const auto parse_endpoint = [&](const std::string& text, std::vector<int>& point) {
			const std::vector<std::string> parts = SplitLegacy(text, ',');
			if (static_cast<int>(parts.size()) != shape.gradients) return false;
			point.clear();
			for (int axis = 0; axis < shape.gradients; ++axis) {
				const std::string token = ToLowerCopy(parts[static_cast<size_t>(axis)]);
				if (token == "var_pos") {
					point.push_back(var_pos);
					continue;
				}
				if (token == "firstlayer") {
					point.push_back(1);
					continue;
				}
				if (token == "lastlayer") {
					point.push_back(shape.layers[axis]);
					continue;
				}
				if (kind == "frozen" && token == "lowerbound") {
					point.push_back(0);
					continue;
				}
				if (kind == "frozen" && token == "upperbound") {
					point.push_back(shape.layers[axis] + 1);
					continue;
				}
				int parsed = std::numeric_limits<int>::min();
				if (!ParseStrict(parts[static_cast<size_t>(axis)], parsed)) return false;
				point.push_back(parsed);
			}
			return true;
		};
		if (!parse_endpoint(endpoints[0], first) || !parse_endpoint(endpoints[1], last)) return false;
	}
	for (int axis = 0; axis < shape.gradients; ++axis) {
		if (first[static_cast<size_t>(axis)] > last[static_cast<size_t>(axis)]) return false;
	}
	coordinates.clear();
	for (int x = first[0]; x <= last[0]; ++x) {
		if (shape.gradients == 1) {
			coordinates.push_back({x});
			continue;
		}
		for (int y = first[1]; y <= last[1]; ++y) {
			if (shape.gradients == 2) {
				coordinates.push_back({x, y});
				continue;
			}
			for (int z = first[2]; z <= last[2]; ++z) coordinates.push_back({x, y, z});
		}
	}
	return true;
}

bool ExtractMaskCoordinates(const json& source,
                            const std::string& key,
                            const ProblemShape& shape,
                            const std::string& base_path,
                            int var_pos,
                            Coordinates& coordinates);

bool ExtractMaskCoordinatesFromJsonFile(const std::string& path,
                                        const std::string& key,
                                        const ProblemShape& shape,
                                        int var_pos,
                                        Coordinates& coordinates) {
	if (!LooksLikeJsonInput(path)) return ParseLegacyMaskFile(path, shape, coordinates);
	json document;
	if (!ReadJsonFile(path, document)) return false;
	if (document.is_object()) {
		if (const auto it = document.find(key); it != document.end()) return ExtractMaskCoordinates(*it, key, shape, path, var_pos, coordinates);
	}
	return ExtractMaskCoordinates(document, key, shape, path, var_pos, coordinates);
}

bool ExtractMaskCoordinates(const json& source,
                            const std::string& key,
                            const ProblemShape& shape,
                            const std::string& base_path,
                            int var_pos,
                            Coordinates& coordinates) {
	if (source.is_string()) {
		const std::string text = source.get<std::string>();
		if (IsLegacyRangeText(text)) return ParseLegacyRange(text, key == "pinned_range" ? "pinned" : "frozen", shape, var_pos, coordinates);
		return ExtractMaskCoordinatesFromJsonFile(ResolveRelativeTo(base_path, text).string(), key, shape, var_pos, coordinates);
	}
	if (ParseCoordinateList(source, shape.gradients, coordinates)) return true;
	if (!source.is_object()) return false;
	if (const json* path = FindMember(source, MASK_FILE_KEYS); path != nullptr && path->is_string()) {
		return ExtractMaskCoordinatesFromJsonFile(ResolveRelativeTo(base_path, path->get<std::string>()).string(), key, shape, var_pos, coordinates);
	}
	if (const auto it = source.find(MASK_COORDINATES_KEY); it != source.end()) return ParseCoordinateList(*it, shape.gradients, coordinates);
	if (const auto it = source.find(key); it != source.end()) return ExtractMaskCoordinates(*it, key, shape, base_path, var_pos, coordinates);
	return false;
}

json CanonicalMask(const Coordinates& coordinates) {
	return {
		{"coordinates", coordinates}
	};
}

bool NormalizeMasks(json& problem, const ProblemShape& shape, const std::string& base_path) {
	const auto mon_section = problem.find("mon");
	if (mon_section == problem.end() || !mon_section->is_object()) return true;
	static constexpr std::array<std::array<const char*, 2>, 2> MASK_KEYS = {{
		{"pinned_range", "pinned_filename"},
		{"frozen_range", "frozen_filename"}
	}};
	for (auto mon = mon_section->begin(); mon != mon_section->end(); ++mon) {
		if (!mon.value().is_object()) continue;
		const std::string freedom = mon.value().value("freedom", std::string{});
		if (freedom == "free") {
			mon.value().erase("pinned_range");
			mon.value().erase("pinned_filename");
			mon.value().erase("frozen_range");
			mon.value().erase("frozen_filename");
		} else if (freedom == "pinned") {
			mon.value().erase("frozen_range");
			mon.value().erase("frozen_filename");
		} else if (freedom == "frozen") {
			mon.value().erase("pinned_range");
			mon.value().erase("pinned_filename");
		}
		const int var_pos = mon.value().value("var_pos", 0);
		const auto normalize = [&](const char* key, const char* filename_key) -> bool {
			if (mon.value().contains(filename_key)) {
				mon.value()[key] = mon.value().at(filename_key);
				mon.value().erase(filename_key);
			}
			const auto it = mon.value().find(key);
			if (it == mon.value().end()) return true;
			Coordinates coordinates;
			if (!ExtractMaskCoordinates(*it, key, shape, base_path, var_pos, coordinates)) return false;
			mon.value()[key] = CanonicalMask(coordinates);
			return true;
		};
		for (const auto& keys : MASK_KEYS) {
			if (!normalize(keys[0], keys[1])) return false;
		}
	}
	return true;
}

bool HasMaskSettings(const json& problem) {
	const auto mon_section = problem.find("mon");
	if (mon_section == problem.end() || !mon_section->is_object()) return false;
	for (auto mon = mon_section->begin(); mon != mon_section->end(); ++mon) {
		if (!mon.value().is_object()) continue;
		if (mon.value().contains("pinned_range") ||
		    mon.value().contains("pinned_filename") ||
		    mon.value().contains("frozen_range") ||
		    mon.value().contains("frozen_filename")) {
			return true;
		}
	}
	return false;
}

bool IsInitialGuessPayload(const json& document) {
	return document.is_object() &&
		document.contains("metadata") &&
		document.contains("monlist") &&
		document.contains("statelist") &&
		document.contains("profiles");
}

const json* FindInitialGuessObject(const json& document) {
	if (IsInitialGuessPayload(document)) return &document;
	if (!document.is_object()) return nullptr;
	if (const auto it = document.find("initial_guess"); it != document.end() && it->is_object()) return &it.value();
	const auto problems = document.find("problems");
	if (problems == document.end()) return nullptr;
	if (problems->is_object()) return FindInitialGuessObject(*problems);
	if (!problems->is_array()) return nullptr;
	for (auto problem = problems->rbegin(); problem != problems->rend(); ++problem) {
		if (!problem->is_object()) continue;
		if (const json* embedded = FindInitialGuessObject(*problem)) return embedded;
	}
	return nullptr;
}

bool LoadInitialGuessSource(const json& source, const std::string& base_path, json& guess) {
	if (source.is_string()) {
		const std::string path = ResolveRelativeTo(base_path, source.get<std::string>()).string();
		json document;
		if (!ReadJsonFile(path, document)) return false;
		const json* embedded = FindInitialGuessObject(document);
		if (embedded == nullptr) return false;
		guess = *embedded;
		return true;
	}
	if (!source.is_object()) return false;
	if (const json* file = FindMember(source, MASK_FILE_KEYS); file != nullptr && file->is_string()) {
		return LoadInitialGuessSource(*file, base_path, guess);
	}
	const json* embedded = FindInitialGuessObject(source);
	if (embedded == nullptr) return false;
	guess = *embedded;
	return true;
}

bool NormalizeInitialGuess(json& problem, const std::string& base_path) {
	auto embed_guess = [&](const json& source) -> bool {
		json guess;
		if (!LoadInitialGuessSource(source, base_path, guess)) return false;
		problem["initial_guess"] = std::move(guess);
		return true;
	};

	if (const auto initial_guess = problem.find("initial_guess"); initial_guess != problem.end()) {
		if (!embed_guess(*initial_guess)) return false;
	}
	const auto sys_section = problem.find("sys");
	if (sys_section != problem.end() && sys_section->is_object()) {
		std::vector<std::string> sys_names;
		for (auto sys = sys_section->begin(); sys != sys_section->end(); ++sys) {
			if (sys.value().is_object()) sys_names.push_back(sys.key());
		}
		for (const std::string& sys_name : sys_names) {
			auto& sys = problem["sys"][sys_name];
			const auto initial_guess = sys.find("initial_guess");
			if (initial_guess == sys.end()) continue;
			if (initial_guess->is_string()) {
				const std::string mode = initial_guess->get<std::string>();
				if (mode == "previous_result" || mode == "file" || mode == "none") continue;
			}
			const json guess_source = *initial_guess;
			if (!embed_guess(guess_source)) return false;
			problem["sys"][sys_name]["initial_guess"] = "file";
		}
	}
	return true;
}

bool NormalizeProblems(json& document, const std::string& base_path) {
	auto normalize_problem = [&](json& problem) {
		if (!problem.is_object()) return true;
		if (HasMaskSettings(problem)) {
			ProblemShape shape;
			if (!LoadProblemShape(problem, shape)) return false;
			if (!NormalizeMasks(problem, shape, base_path)) return false;
		}
		return NormalizeInitialGuess(problem, base_path);
	};
	if (document.is_array()) {
		for (auto& problem : document) {
			if (!normalize_problem(problem)) return false;
		}
		return true;
	}
	if (!document.is_object()) return true;
	const auto problems = document.find("problems");
	if (problems == document.end()) return normalize_problem(document);
	if (!problems->is_array()) return normalize_problem(document["problems"]);
	for (auto& problem : document["problems"]) {
		if (!normalize_problem(problem)) return false;
	}
	return true;
}

} // namespace

namespace io::input {

bool PrepareInputFile(const std::string& requested_path, std::string& json_path) {
	json document;
	if (LooksLikeJsonInput(requested_path)) {
		if (!ReadJsonFile(requested_path, document)) return false;
	} else {
		std::vector<std::string> entries;
		if (!ReadLegacyEntries(requested_path, entries)) return false;
		json problems = json::array();
		if (!ConvertLegacyEntries(entries, problems)) return false;
		document = {
			{"problems", std::move(problems)}
		};
	}
	if (!NormalizeProblems(document, requested_path)) return false;
	const std::filesystem::path emitted = EmittedJsonPath(requested_path);
	std::ofstream out(emitted.c_str(), std::ios::out | std::ios::trunc);
	if (!out.is_open()) {
		std::cout << "Failed to write preprocessed input file " << emitted.string() << std::endl;
		return false;
	}
	out << document.dump(2) << std::endl;
	json_path = emitted.string();
	return true;
}

} // namespace io::input
