#ifndef INPUTxH
#define INPUTxH
#include "namics.h"
#include <cctype>

inline std::string ToLowerCopy(std::string value) {
	std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
		return static_cast<char>(std::tolower(c));
	});
	return value;
}

template <typename T>
inline bool ParseStrict(const std::string& s, T& value) {
	if (s.empty()) return false;
	std::istringstream stream(s);
	stream >> std::noskipws >> value;
	return !stream.fail() && stream.eof();
}

template <>
inline bool ParseStrict<bool>(const std::string& s, bool& value) {
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

template <typename T>
inline bool ParseOrReport(const std::string& s, T& value, const std::string& error) {
	if (ParseStrict(s, value)) return true;
	std::cout << error << std::endl;
	return false;
}

template <typename T, typename L, typename H>
inline bool ParseRangeOrReport(const std::string& s, T& value, L low, H high, const std::string& error) {
	if (!ParseStrict(s, value)) {
		std::cout << error << std::endl;
		return false;
	}
	const T low_value = static_cast<T>(low);
	const T high_value = static_cast<T>(high);
	if (value < low_value || value > high_value) {
		std::cout << "Value out of range: " << error << std::endl;
		return false;
	}
	return true;
}

template <typename T>
inline bool ContainsValue(const std::vector<T>& values, const T& target, int* pos = nullptr) {
	const auto it = std::find(values.begin(), values.end(), target);
	if (it == values.end()) return false;
	if (pos != nullptr) *pos = static_cast<int>(std::distance(values.begin(), it));
	return true;
}

inline int ParseInt(const std::string& s, int fallback) {
	int value = fallback;
	ParseStrict(s, value);
	return value;
}

inline bool ParseInt(const std::string& s, int& value, const std::string& error) {
	return ParseOrReport(s, value, error);
}

template <typename L, typename H>
inline bool ParseInt(const std::string& s, int& value, L low, H high, const std::string& error) {
	return ParseRangeOrReport(s, value, low, high, error);
}

inline Real ParseReal(const std::string& s, Real fallback) {
	Real value = fallback;
	ParseStrict(s, value);
	return value;
}

inline bool ParseReal(const std::string& s, Real& value, const std::string& error) {
	return ParseOrReport(s, value, error);
}

template <typename L, typename H>
inline bool ParseReal(const std::string& s, Real& value, L low, H high, const std::string& error) {
	return ParseRangeOrReport(s, value, low, high, error);
}

inline bool ParseBool(const std::string& s, bool fallback) {
	bool value = fallback;
	ParseStrict(s, value);
	return value;
}

inline bool ParseBool(const std::string& s, bool& value, const std::string& error) {
	return ParseOrReport(s, value, error);
}

inline std::string ParseString(const std::string& s, const std::string& fallback) {
	return s.empty() ? fallback : s;
}

inline bool ParseString(const std::string& s, std::string& value, const std::string& error) {
	if (s.empty()) {
		std::cout << error << std::endl;
		return false;
	}
	value = s;
	return true;
}

inline bool ParseString(const std::string& s, std::string& value, const std::vector<std::string>& allowed, const std::string& error) {
	if (!ParseString(s, value, error)) return false;
	if (!ContainsValue(allowed, value)) {
		std::cout << error << " value '" << value << "' is not allowed. Select from: " << std::endl;
		for (const std::string& item : allowed) std::cout << item << " ; ";
		std::cout << std::endl;
		return false;
	}
	return true;
}

class Input {
public:
	Input(const std::string&);

~Input();

	std::string name;
	std::ifstream in_file;
	std::ifstream inc_file;


	std::string In_buffer;
	std::string string_value;
	bool Input_error;
	std::string filename;
	std::string output_path;

	std::vector<std::string> KEYS;
	std::vector<std::string> SysList;
	std::vector<std::string> MolList;
	std::vector<std::string> MonList;
	std::vector<std::string> LatList;
	std::vector<std::string> NewtonList;
	std::vector<std::string> OutputList;
	std::vector<std::string> elems;
	std::vector<std::string> StateList;
	std::vector<std::string> ReactionList;


	void PrintList(const std::vector<std::string>&) const;
	std::vector<std::string>& split(const std::string&, char, std::vector<std::string>&) const;
	bool TestNum(std::vector<std::string>&, const std::string&, int, int, int) const;
	int GetNumStarts(void) const;
	bool CheckParameters(const std::string&, const std::string&, int, const std::vector<std::string>&, ParameterStore&) const;
	bool LoadItems(const std::string&, std::vector<std::string>&, std::vector<std::string>&, std::vector<std::string>&) const;
	bool CheckInput(void);
	bool InSet(const std::vector<std::string>&, const std::string&) const;
	bool InSet(const std::vector<int>&, int) const;
	bool InSet(const std::vector<int>&, int&, int) const;
	bool InSet(const std::vector<std::string>&, int&, const std::string&) const;
	bool ReadFile(const std::string&, std::string&) const;
	bool EvenBrackets(const std::string&, std::vector<int>&, std::vector<int>&) const;
	bool EvenSquareBrackets(const std::string&, std::vector<int>&, std::vector<int>&) const;
	bool MakeLists(int);
	const std::string& GetOutputPath() const;
	std::string ResolvePath(const std::string&) const;

private:
	void parseOutputInfo();
	bool OutputPathExists() const;
};

#endif
