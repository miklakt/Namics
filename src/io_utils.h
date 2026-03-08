#ifndef IO_UTILSxH
#define IO_UTILSxH

#include "namics.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace io {

namespace detail {

inline int GuessProfileSize(int mx, int my, int mz, int fjc) {
	if (my == 0) return (mx + 2 * fjc);
	if (mz == 0) return (mx + 2 * fjc) * (my + 2 * fjc);
	return (mx + 2 * fjc) * (my + 2 * fjc) * (mz + 2 * fjc);
}

inline bool ReadWholeFile(const std::string& filename, std::string& content) {
	std::ifstream input(filename.c_str());
	if (!input.is_open()) return false;
	std::ostringstream buffer;
	buffer << input.rdbuf();
	content = buffer.str();
	return static_cast<bool>(input) || input.eof();
}

inline std::string JsonUnescape(const std::string& value) {
	std::string out;
	out.reserve(value.size());
	bool escaped = false;
	for (const char c : value) {
		if (escaped) {
			switch (c) {
				case 'n': out.push_back('\n'); break;
				case 'r': out.push_back('\r'); break;
				case 't': out.push_back('\t'); break;
				case '\\': out.push_back('\\'); break;
				case '"': out.push_back('"'); break;
				default: out.push_back(c); break;
			}
			escaped = false;
		} else if (c == '\\') {
			escaped = true;
		} else {
			out.push_back(c);
		}
	}
	if (escaped) out.push_back('\\');
	return out;
}

inline size_t SkipWs(const std::string& text, size_t pos) {
	while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
	return pos;
}

inline size_t FindMatchingDelimiter(const std::string& text, size_t open_pos, char open_ch, char close_ch) {
	if (open_pos >= text.size() || text[open_pos] != open_ch) return std::string::npos;
	int depth = 0;
	bool in_string = false;
	bool escaped = false;
	for (size_t i = open_pos; i < text.size(); ++i) {
		const char c = text[i];
		if (in_string) {
			if (escaped) {
				escaped = false;
			} else if (c == '\\') {
				escaped = true;
			} else if (c == '"') {
				in_string = false;
			}
			continue;
		}
		if (c == '"') {
			in_string = true;
			continue;
		}
		if (c == open_ch) {
			++depth;
			continue;
		}
		if (c == close_ch) {
			--depth;
			if (depth == 0) return i;
		}
	}
	return std::string::npos;
}

inline bool ExtractEnclosedForKey(const std::string& object_text,
                                  const std::string& key,
                                  char open_ch,
                                  char close_ch,
                                  std::string& enclosed) {
	const std::string needle = "\"" + key + "\"";
	const size_t key_pos = object_text.find(needle);
	if (key_pos == std::string::npos) return false;
	const size_t colon = object_text.find(':', key_pos + needle.size());
	if (colon == std::string::npos) return false;
	const size_t value_pos = SkipWs(object_text, colon + 1);
	if (value_pos >= object_text.size() || object_text[value_pos] != open_ch) return false;
	const size_t end_pos = FindMatchingDelimiter(object_text, value_pos, open_ch, close_ch);
	if (end_pos == std::string::npos) return false;
	enclosed = object_text.substr(value_pos, end_pos - value_pos + 1);
	return true;
}

inline bool ExtractStringForKey(const std::string& object_text,
                                const std::string& key,
                                std::string& value) {
	const std::string needle = "\"" + key + "\"";
	const size_t key_pos = object_text.find(needle);
	if (key_pos == std::string::npos) return false;
	const size_t colon = object_text.find(':', key_pos + needle.size());
	if (colon == std::string::npos) return false;
	size_t pos = SkipWs(object_text, colon + 1);
	if (pos >= object_text.size() || object_text[pos] != '"') return false;
	++pos;
	std::string encoded;
	bool escaped = false;
	while (pos < object_text.size()) {
		const char c = object_text[pos++];
		if (escaped) {
			encoded.push_back('\\');
			encoded.push_back(c);
			escaped = false;
			continue;
		}
		if (c == '\\') {
			escaped = true;
			continue;
		}
		if (c == '"') {
			value = JsonUnescape(encoded);
			return true;
		}
		encoded.push_back(c);
	}
	return false;
}

inline bool ExtractTokenForKey(const std::string& object_text,
                               const std::string& key,
                               std::string& token) {
	const std::string needle = "\"" + key + "\"";
	const size_t key_pos = object_text.find(needle);
	if (key_pos == std::string::npos) return false;
	const size_t colon = object_text.find(':', key_pos + needle.size());
	if (colon == std::string::npos) return false;
	const size_t start = SkipWs(object_text, colon + 1);
	if (start >= object_text.size()) return false;
	size_t end = start;
	while (end < object_text.size()) {
		const char c = object_text[end];
		if (c == ',' || c == '}' || c == '\n' || c == '\r' || c == '\t' || c == ' ') break;
		++end;
	}
	if (end <= start) return false;
	token = object_text.substr(start, end - start);
	return true;
}

inline bool ExtractIntForKey(const std::string& object_text, const std::string& key, int& value) {
	std::string token;
	if (!ExtractTokenForKey(object_text, key, token)) return false;
	char* end_ptr = nullptr;
	const long parsed = std::strtol(token.c_str(), &end_ptr, 10);
	if (end_ptr == token.c_str() || *end_ptr != '\0') return false;
	value = static_cast<int>(parsed);
	return true;
}

inline bool ExtractBoolForKey(const std::string& object_text, const std::string& key, bool& value) {
	std::string token;
	if (!ExtractTokenForKey(object_text, key, token)) return false;
	if (token == "true") {
		value = true;
		return true;
	}
	if (token == "false") {
		value = false;
		return true;
	}
	return false;
}

inline bool ParseStringArray(const std::string& array_text, std::vector<std::string>& values) {
	values.clear();
	if (array_text.size() < 2 || array_text.front() != '[' || array_text.back() != ']') return false;
	size_t pos = 1;
	while (pos < array_text.size() - 1) {
		pos = SkipWs(array_text, pos);
		if (pos >= array_text.size() - 1) break;
		if (array_text[pos] == ',') {
			++pos;
			continue;
		}
		if (array_text[pos] != '"') return false;
		++pos;
		std::string encoded;
		bool escaped = false;
		while (pos < array_text.size() - 1) {
			const char c = array_text[pos++];
			if (escaped) {
				encoded.push_back('\\');
				encoded.push_back(c);
				escaped = false;
				continue;
			}
			if (c == '\\') {
				escaped = true;
				continue;
			}
			if (c == '"') {
				values.push_back(JsonUnescape(encoded));
				break;
			}
			encoded.push_back(c);
		}
	}
	return true;
}

template <typename RealT>
inline bool ParseNumberArray(const std::string& array_text, std::vector<RealT>& values) {
	values.clear();
	if (array_text.size() < 2 || array_text.front() != '[' || array_text.back() != ']') return false;
	const char* ptr = array_text.c_str() + 1;
	const char* end = array_text.c_str() + array_text.size() - 1;
	while (ptr < end) {
		while (ptr < end && (std::isspace(static_cast<unsigned char>(*ptr)) || *ptr == ',')) ++ptr;
		if (ptr >= end) break;
		char* next = nullptr;
		const long double parsed = std::strtold(ptr, &next);
		if (next == ptr) return false;
		values.push_back(static_cast<RealT>(parsed));
		ptr = next;
	}
	return true;
}

template <typename RealT>
inline bool ExtractNumberArrayForKey(const std::string& object_text,
                                     const std::string& key,
                                     std::vector<RealT>& values) {
	std::string array_text;
	if (!ExtractEnclosedForKey(object_text, key, '[', ']', array_text)) return false;
	return ParseNumberArray(array_text, values);
}

inline bool ExtractLastEmbeddedInitialGuessObject(const std::string& json_text, std::string& guess_object) {
	const std::string key = "\"initial_guess\"";
	size_t key_pos = json_text.rfind(key);
	while (key_pos != std::string::npos) {
		const size_t colon = json_text.find(':', key_pos + key.size());
		if (colon != std::string::npos) {
			const size_t value_pos = SkipWs(json_text, colon + 1);
			if (value_pos < json_text.size() && json_text[value_pos] == '{') {
				const size_t end_pos = FindMatchingDelimiter(json_text, value_pos, '{', '}');
				if (end_pos != std::string::npos) {
					guess_object = json_text.substr(value_pos, end_pos - value_pos + 1);
					return true;
				}
			}
		}
		if (key_pos == 0) break;
		key_pos = json_text.rfind(key, key_pos - 1);
	}
	return false;
}

template <typename RealT>
inline bool ParseInitialGuessJsonObject(const std::string& object_text,
                                        RealT* x,
                                        std::string& method,
                                        std::vector<std::string>& monlist,
                                        std::vector<std::string>& statelist,
                                        bool& charged,
                                        int& mx,
                                        int& my,
                                        int& mz,
                                        int& fjc,
                                        int readx) {
	std::string metadata_object;
	if (!ExtractEnclosedForKey(object_text, "metadata", '{', '}', metadata_object)) return false;
	if (!ExtractStringForKey(metadata_object, "method", method)) return false;
	if (!ExtractIntForKey(metadata_object, "mx", mx)) return false;
	if (!ExtractIntForKey(metadata_object, "my", my)) return false;
	if (!ExtractIntForKey(metadata_object, "mz", mz)) return false;
	if (!ExtractIntForKey(metadata_object, "fjc", fjc)) return false;
	if (!ExtractBoolForKey(metadata_object, "charged", charged)) return false;

	std::string monlist_array;
	if (!ExtractEnclosedForKey(object_text, "monlist", '[', ']', monlist_array)) return false;
	if (!ParseStringArray(monlist_array, monlist)) return false;

	std::string statelist_array;
	if (!ExtractEnclosedForKey(object_text, "statelist", '[', ']', statelist_array)) return false;
	if (!ParseStringArray(statelist_array, statelist)) return false;

	if (readx == 0) return true;

	std::string profiles_object;
	if (!ExtractEnclosedForKey(object_text, "profiles", '{', '}', profiles_object)) return false;

	const int m = GuessProfileSize(mx, my, mz, fjc);
	const int profile_count = static_cast<int>(monlist.size() + statelist.size() + (charged ? 1 : 0));
	std::vector<RealT> assembled;
	assembled.reserve(static_cast<size_t>(profile_count) * static_cast<size_t>(m));

	std::vector<RealT> values;
	for (const std::string& mon_name : monlist) {
		if (!ExtractNumberArrayForKey(profiles_object, "mon:" + mon_name, values)) return false;
		if (static_cast<int>(values.size()) != m) return false;
		assembled.insert(assembled.end(), values.begin(), values.end());
	}
	for (const std::string& state_name : statelist) {
		if (!ExtractNumberArrayForKey(profiles_object, "state:" + state_name, values)) return false;
		if (static_cast<int>(values.size()) != m) return false;
		assembled.insert(assembled.end(), values.begin(), values.end());
	}
	if (charged) {
		if (!ExtractNumberArrayForKey(profiles_object, "psi", values)) return false;
		if (static_cast<int>(values.size()) != m) return false;
		assembled.insert(assembled.end(), values.begin(), values.end());
	}

	for (size_t i = 0; i < assembled.size(); ++i) x[i] = assembled[i];
	return true;
}

template <typename RealT>
inline bool ReadInitialGuessJson(const std::string& filename,
                                 RealT* x,
                                 std::string& method,
                                 std::vector<std::string>& monlist,
                                 std::vector<std::string>& statelist,
                                 bool& charged,
                                 int& mx,
                                 int& my,
                                 int& mz,
                                 int& fjc,
                                 int readx) {
	std::string content;
	if (!ReadWholeFile(filename, content)) {
		std::cout << "inputfile " << filename << " is not found. Read guess for initial guess failed" << std::endl;
		return false;
	}
	if (ParseInitialGuessJsonObject(content, x, method, monlist, statelist, charged, mx, my, mz, fjc, readx)) {
		return true;
	}

	std::string embedded_guess;
	if (!ExtractLastEmbeddedInitialGuessObject(content, embedded_guess)) {
		std::cout << "No initial_guess found in " << filename << ". Read guess for initial guess failed" << std::endl;
		return false;
	}
	if (!ParseInitialGuessJsonObject(embedded_guess, x, method, monlist, statelist, charged, mx, my, mz, fjc, readx)) {
		std::cout << "Found initial_guess in " << filename << " but failed to parse it. Read guess for initial guess failed" << std::endl;
		return false;
	}
	return true;
}

} // namespace detail

inline bool ReadSanitizedFile(const std::string& filename, std::string& buffer) {
	std::ifstream input(filename.c_str());
	if (!input.is_open()) {
		std::cout << "Inputfile " << filename << " is not found. " << std::endl;
		return false;
	}

	buffer.clear();
	std::string line;
	while (std::getline(input, line)) {
		line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
		line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
		if (line.empty()) continue;
		if (line.size() >= 2 && line.substr(0, 2) == "//") continue;
		buffer.append(line).append("#");
	}

	if (buffer.empty()) {
		std::cout << "File " << filename << " is empty " << std::endl;
		return false;
	}
	return true;
}

template <typename RealT>
inline bool ReadExternalPotentialJson(const std::string& filename,
                                      std::vector<RealT>& values) {
	std::string content;
	if (!detail::ReadWholeFile(filename, content)) {
		std::cout << "Inputfile " << filename << " is not found. " << std::endl;
		return false;
	}

	const auto try_object = [&](const std::string& object_text) {
		if (detail::ExtractNumberArrayForKey(object_text, "external_potential", values)) return true;
		std::string profiles_object;
		if (detail::ExtractEnclosedForKey(object_text, "profiles", '{', '}', profiles_object) &&
		    detail::ExtractNumberArrayForKey(profiles_object, "external_potential", values)) {
			return true;
		}
		return false;
	};

	const std::string key = "\"problems\"";
	const size_t key_pos = content.find(key);
	if (key_pos != std::string::npos) {
		const size_t colon = content.find(':', key_pos + key.size());
		if (colon != std::string::npos) {
			const size_t value_pos = detail::SkipWs(content, colon + 1);
			if (value_pos < content.size() && content[value_pos] == '[') {
				const size_t array_end = detail::FindMatchingDelimiter(content, value_pos, '[', ']');
				if (array_end != std::string::npos) {
					size_t pos = value_pos + 1;
					std::string object_text;
					while (pos < array_end) {
						pos = detail::SkipWs(content, pos);
						if (pos >= array_end) break;
						if (content[pos] == ',') {
							++pos;
							continue;
						}
						if (content[pos] != '{') break;
						const size_t object_end = detail::FindMatchingDelimiter(content, pos, '{', '}');
						if (object_end == std::string::npos || object_end > array_end) break;
						object_text = content.substr(pos, object_end - pos + 1);
						pos = object_end + 1;
					}
					if (!object_text.empty() && try_object(object_text)) return true;
				}
			}
		}
	}

	if (try_object(content)) return true;

	std::cout << "Unable to find json array 'external_potential' in " << filename << std::endl;
	return false;
}

template <typename RealT>
inline bool ReadInitialGuess(const std::string& filename,
                             RealT* x,
                             std::string& method,
                             std::vector<std::string>& monlist,
                             std::vector<std::string>& statelist,
                             bool& charged,
                             int& mx,
                             int& my,
                             int& mz,
                             int& fjc,
                             int readx) {
	return detail::ReadInitialGuessJson(filename, x, method, monlist, statelist, charged, mx, my, mz, fjc, readx);
}

} // namespace io

#endif
