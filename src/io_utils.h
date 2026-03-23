#ifndef IO_UTILSxH
#define IO_UTILSxH

#include <algorithm>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <span>
#include <string>
#include <vector>

namespace io {

namespace detail {

using json = nlohmann::ordered_json;

inline bool ReadJsonFile(const std::string& filename, json& document) {
	std::ifstream input(filename.c_str());
	if (!input.is_open()) {
		std::cout << "Inputfile " << filename << " is not found. " << std::endl;
		return false;
	}
	try {
		document = json::parse(input, nullptr, true, true);
	} catch (const std::exception& error) {
		std::cout << "Failed to parse JSON file " << filename << ": " << error.what() << std::endl;
		return false;
	}
	return true;
}

inline const json* FindLastProblemObject(const json& document) {
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
	if (!document.is_object()) return nullptr;
	if (const auto guess = document.find("initial_guess"); guess != document.end() && guess->is_object()) return &*guess;
	if (const auto profiles = document.find("profiles"); profiles != document.end() && profiles->is_object()) return &document;
	if (const auto* problem = FindLastProblemObject(document)) return FindInitialGuessObject(*problem);
	return nullptr;
}

inline const json* FindExternalPotentialNode(const json& document) {
	if (document.is_array()) return &document;
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

template <typename RealT>
inline bool ReadInitialGuess(const std::string& filename,
                             std::span<RealT> x,
                             const std::vector<std::string>& monlist,
                             const std::vector<std::string>& statelist,
                             bool charged) {
	detail::json document;
	if (!detail::ReadJsonFile(filename, document)) {
		std::cout << "Read guess for initial guess failed" << std::endl;
		return false;
	}
	const auto* guess = detail::FindInitialGuessObject(document);
	if (guess == nullptr) {
		std::cout << "No initial_guess found in " << filename << ". Read guess for initial guess failed" << std::endl;
		return false;
	}
	const int count = static_cast<int>(monlist.size() + statelist.size() + (charged ? 1 : 0));
	if (count == 0) return x.empty();
	if (static_cast<int>(x.size()) % count != 0) return false;
	return detail::ReadInitialGuessProfiles(*guess, x, monlist, statelist, charged, static_cast<int>(x.size()) / count);
}

} // namespace io

#endif
