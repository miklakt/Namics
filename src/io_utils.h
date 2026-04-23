#ifndef IO_UTILSxH
#define IO_UTILSxH

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <span>
#include <string>
#include <vector>

namespace io {

namespace detail {

using json = nlohmann::ordered_json;

inline void SplitJsonSubpath(const std::string& path, std::string& filename, std::string& pointer) {
	const std::string marker = ".json/";
	const size_t pos = path.rfind(marker);
	if (pos == std::string::npos) {
		filename = path;
		pointer.clear();
		return;
	}
	filename = path.substr(0, pos + 5);
	pointer = path.substr(pos + 5);
}

inline bool ReadJsonFile(const std::string& filename, json& document) {
	std::string json_filename;
	std::string json_pointer;
	SplitJsonSubpath(filename, json_filename, json_pointer);
	std::ifstream input(json_filename.c_str());
	if (!input.is_open()) {
		std::cout << "Inputfile " << json_filename << " is not found. " << std::endl;
		return false;
	}
	try {
		document = json::parse(input, nullptr, true, true);
	} catch (const std::exception& error) {
		std::cout << "Failed to parse JSON file " << json_filename << ": " << error.what() << std::endl;
		return false;
	}
	if (!json_pointer.empty()) {
		try {
			document = document.at(json::json_pointer(json_pointer));
		} catch (const std::exception& error) {
			std::cout << "Failed to select JSON subtree '" << json_pointer << "' in " << json_filename << ": " << error.what() << std::endl;
			return false;
		}
	}
	return true;
}

inline const json* FindLastProblemObject(const json& document) {
	if (document.is_array()) {
		for (auto it = document.rbegin(); it != document.rend(); ++it) {
			if (it->is_object()) return &*it;
		}
		return nullptr;
	}
	if (!document.is_object()) return nullptr;
	const auto problems = document.find("problems");
	if (problems == document.end()) return nullptr;
	if (problems->is_object()) return &*problems;
	if (!problems->is_array()) return nullptr;
	for (auto it = problems->rbegin(); it != problems->rend(); ++it) {
		if (it->is_object()) return &*it;
	}
	return nullptr;
}

inline const json* FindInitialGuessObject(const json& document) {
	if (document.is_array()) {
		for (auto it = document.rbegin(); it != document.rend(); ++it) {
			if (const auto* guess = FindInitialGuessObject(*it)) return guess;
		}
		return nullptr;
	}
	if (!document.is_object()) return nullptr;
	if (const auto guess = document.find("initial_guess"); guess != document.end() && guess->is_object()) return &*guess;
	if (const auto profiles = document.find("profiles"); profiles != document.end() && profiles->is_object()) return &document;
	if (const auto* problem = FindLastProblemObject(document)) return FindInitialGuessObject(*problem);
	return nullptr;
}

inline const json* FindExternalPotentialNode(const json& document) {
	if (document.is_array()) {
		if (std::none_of(document.begin(), document.end(), [](const json& item) { return item.is_object(); })) return &document;
		for (auto it = document.rbegin(); it != document.rend(); ++it) {
			if (const auto* node = FindExternalPotentialNode(*it)) return node;
		}
		return nullptr;
	}
	if (document.is_object()) {
		if (const auto values = document.find("external_potential"); values != document.end() && values->is_array()) return &*values;
		if (const auto profiles = document.find("profiles"); profiles != document.end() && profiles->is_object()) {
			const auto values = profiles->find("external_potential");
			if (values != profiles->end() && values->is_array()) return &*values;
		}
		if (const auto* problem = FindLastProblemObject(document)) return FindExternalPotentialNode(*problem);
	}
	return nullptr;
}

template <typename RealT>
inline bool ReadDenseArray(const json& node, std::vector<RealT>& values) {
	if (!node.is_array()) return false;
	try {
		values = node.template get<std::vector<RealT>>();
	} catch (const json::exception&) {
		values.clear();
		return false;
	}
	return true;
}

template <typename RealT>
inline bool ReadProfileArray(const json& node, int profile_size, std::vector<RealT>& values) {
	if (ReadDenseArray(node, values)) return static_cast<int>(values.size()) == profile_size;
	if (!node.is_object()) return false;
	if (const auto dense = node.find("dense"); dense != node.end()) return ReadProfileArray(*dense, profile_size, values);
	if (const auto dense = node.find("values"); dense != node.end()) return ReadProfileArray(*dense, profile_size, values);
	const auto sparse = node.find("sparse");
	if (sparse == node.end() || !sparse->is_array()) return false;
	values.assign(static_cast<size_t>(profile_size), 0);
	for (const auto& entry : *sparse) {
		int index = -1;
		RealT value = 0;
		try {
			if (entry.is_array() && entry.size() == 2) {
				index = entry[0].get<int>();
				value = entry[1].get<RealT>();
			} else if (entry.is_object()) {
				index = entry.at("index").get<int>();
				value = entry.at("value").get<RealT>();
			} else {
				return false;
			}
		} catch (const json::exception&) {
			return false;
		}
		if (index < 0 || index >= profile_size) return false;
		values[static_cast<size_t>(index)] = value;
	}
	return true;
}

template <typename RealT>
inline bool FlattenProfileArray(const json& node, std::vector<RealT>& values) {
	if (node.is_array()) {
		for (const auto& item : node) {
			if (!FlattenProfileArray(item, values)) return false;
		}
		return true;
	}
	try {
		values.push_back(node.template get<RealT>());
		return true;
	} catch (const json::exception&) {
		return false;
	}
}

inline const json* FindProfileInOutputSection(const json& section, const std::string& name, const std::string& prop) {
	const auto item = section.find(name);
	if (item == section.end() || !item->is_object()) return nullptr;
	const auto profile = item->find(prop);
	return profile != item->end() ? &*profile : nullptr;
}

template <typename RealT, typename PadProfileFn>
inline bool ReadOutputProfileArray(const json& node, int profile_size, std::vector<RealT>& values, PadProfileFn&& pad_profile) {
	values.clear();
	if (!FlattenProfileArray(node, values)) return false;
	return static_cast<int>(values.size()) == profile_size || pad_profile(values);
}

template <typename RealT>
inline bool WriteInitialGuessProfiles(json& guess,
                                     std::span<const RealT> values,
                                     const std::vector<std::string>& monlist,
                                     const std::vector<std::string>& statelist,
                                     bool charged,
                                     int profile_size) {
	const int count = static_cast<int>(monlist.size() + statelist.size() + (charged ? 1 : 0));
	if (static_cast<int>(values.size()) != count * profile_size) return false;
	guess = json::object();
	guess["profiles"] = json::object();
	size_t offset = 0;
	const auto write_profile = [&](const std::string& key) {
		guess["profiles"][key] = std::vector<RealT>(
			values.begin() + static_cast<std::ptrdiff_t>(offset * static_cast<size_t>(profile_size)),
			values.begin() + static_cast<std::ptrdiff_t>((offset + 1) * static_cast<size_t>(profile_size)));
		++offset;
	};
	for (const std::string& mon_name : monlist) write_profile("mon:" + mon_name);
	for (const std::string& state_name : statelist) write_profile("state:" + state_name);
	if (charged) write_profile("psi");
	return true;
}

template <typename RealT>
inline bool ReadInitialGuessProfiles(const json& guess,
                                     std::span<RealT> x,
                                     const std::vector<std::string>& monlist,
                                     const std::vector<std::string>& statelist,
                                     bool charged,
                                     int profile_size) {
	const auto profiles = guess.find("profiles");
	if (profiles == guess.end() || !profiles->is_object()) return false;
	const int count = static_cast<int>(monlist.size() + statelist.size() + (charged ? 1 : 0));
	if (static_cast<int>(x.size()) != count * profile_size) return false;
	std::vector<RealT> values;
	size_t offset = 0;
	for (const std::string& mon_name : monlist) {
		const auto profile = profiles->find("mon:" + mon_name);
		if (profile == profiles->end() || !ReadProfileArray(*profile, profile_size, values)) return false;
		std::copy(values.begin(), values.end(), x.begin() + static_cast<std::ptrdiff_t>(offset));
		offset += static_cast<size_t>(profile_size);
	}
	for (const std::string& state_name : statelist) {
		const auto profile = profiles->find("state:" + state_name);
		if (profile == profiles->end() || !ReadProfileArray(*profile, profile_size, values)) return false;
		std::copy(values.begin(), values.end(), x.begin() + static_cast<std::ptrdiff_t>(offset));
		offset += static_cast<size_t>(profile_size);
	}
	if (charged) {
		const auto profile = profiles->find("psi");
		if (profile == profiles->end() || !ReadProfileArray(*profile, profile_size, values)) return false;
		std::copy(values.begin(), values.end(), x.begin() + static_cast<std::ptrdiff_t>(offset));
	}
	return true;
}

} // namespace detail

template <typename RealT>
inline bool ReadExternalPotentialJson(const std::string& filename, std::vector<RealT>& values) {
	detail::json document;
	if (!detail::ReadJsonFile(filename, document)) return false;
	const auto* node = detail::FindExternalPotentialNode(document);
	if (node != nullptr && detail::ReadDenseArray(*node, values)) return true;
	std::cout << "Unable to find json array 'external_potential' in " << filename << std::endl;
	return false;
}

template <typename RealT, typename ResolveProfileFn>
inline bool ReadInitialGuess(const std::string& filename,
                             std::span<RealT> x,
                             const std::vector<std::string>& monlist,
                             const std::vector<std::string>& statelist,
                             bool charged,
                             ResolveProfileFn&& resolve_profile) {
	const int count = static_cast<int>(monlist.size() + statelist.size() + (charged ? 1 : 0));
	detail::json document;
	if (!detail::ReadJsonFile(filename, document)) {
		std::cout << "Read guess for initial guess failed" << std::endl;
		return false;
	}
	if (const auto* guess = detail::FindInitialGuessObject(document); guess != nullptr) {
		if (count == 0) return x.empty();
		if (static_cast<int>(x.size()) % count != 0) return false;
		return detail::ReadInitialGuessProfiles(*guess, x, monlist, statelist, charged, static_cast<int>(x.size()) / count);
	}
	if (count == 0) return x.empty();

	detail::json guess;
	guess["profiles"] = detail::json::object();
	std::vector<RealT> values;
	for (const std::string& mon_name : monlist) {
		if (!resolve_profile(document, "mon:" + mon_name, values)) {
			std::cout << "No initial_guess found in " << filename << ". Read guess for initial guess failed" << std::endl;
			return false;
		}
		guess["profiles"]["mon:" + mon_name] = values;
	}
	for (const std::string& state_name : statelist) {
		if (!resolve_profile(document, "state:" + state_name, values)) {
			std::cout << "No initial_guess found in " << filename << ". Read guess for initial guess failed" << std::endl;
			return false;
		}
		guess["profiles"]["state:" + state_name] = values;
	}
	if (charged) {
		if (!resolve_profile(document, "psi", values)) {
			std::cout << "No initial_guess found in " << filename << ". Read guess for initial guess failed" << std::endl;
			return false;
		}
		guess["profiles"]["psi"] = values;
	}
	return detail::ReadInitialGuessProfiles(guess, x, monlist, statelist, charged, static_cast<int>(x.size()) / count);
}

} // namespace io

#endif
