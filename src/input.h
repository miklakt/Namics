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
inline bool ParseStrict(const string& s, T& value) {
	if (s.empty()) return false;
	std::istringstream stream(s);
	stream >> std::noskipws >> value;
	return !stream.fail() && stream.eof();
}

template <>
inline bool ParseStrict<bool>(const string& s, bool& value) {
	const string lowered = ToLowerCopy(s);
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
inline bool ParseOrReport(const string& s, T& value, const std::string& error) {
	if (ParseStrict(s, value)) return true;
	cout << error << endl;
	return false;
}

template <typename T, typename L, typename H>
inline bool ParseRangeOrReport(const string& s, T& value, L low, H high, const std::string& error) {
	if (!ParseStrict(s, value)) {
		cout << error << endl;
		return false;
	}
	const T low_value = static_cast<T>(low);
	const T high_value = static_cast<T>(high);
	if (value < low_value || value > high_value) {
		cout << "Value out of range: " << error << endl;
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

inline int ParseInt(const string& s, int fallback) {
	int value = fallback;
	ParseStrict(s, value);
	return value;
}

inline bool ParseInt(const string& s, int& value, const std::string& error) {
	return ParseOrReport(s, value, error);
}

template <typename L, typename H>
inline bool ParseInt(const string& s, int& value, L low, H high, const std::string& error) {
	return ParseRangeOrReport(s, value, low, high, error);
}

inline Real ParseReal(const string& s, Real fallback) {
	Real value = fallback;
	ParseStrict(s, value);
	return value;
}

inline bool ParseReal(const string& s, Real& value, const std::string& error) {
	return ParseOrReport(s, value, error);
}

template <typename L, typename H>
inline bool ParseReal(const string& s, Real& value, L low, H high, const std::string& error) {
	return ParseRangeOrReport(s, value, low, high, error);
}

inline bool ParseBool(const string& s, bool fallback) {
	bool value = fallback;
	ParseStrict(s, value);
	return value;
}

inline bool ParseBool(const string& s, bool& value, const std::string& error) {
	return ParseOrReport(s, value, error);
}

inline string ParseString(const string& s, const string& fallback) {
	return s.empty() ? fallback : s;
}

inline bool ParseString(const string& s, string& value, const std::string& error) {
	if (s.empty()) {
		cout << error << endl;
		return false;
	}
	value = s;
	return true;
}

inline bool ParseString(const string& s, string& value, const std::vector<std::string>& allowed, const std::string& error) {
	if (!ParseString(s, value, error)) return false;
	if (!ContainsValue(allowed, value)) {
		cout << error << " value '" << value << "' is not allowed. Select from: " << endl;
		for (const std::string& item : allowed) cout << item << " ; ";
		cout << endl;
		return false;
	}
	return true;
}

class Input {
public:
	Input(const string&);

~Input();

	string name;
	ifstream in_file;
	ifstream inc_file;


	std::string In_buffer;
	string string_value;
	bool Input_error;
	string filename;
	std::string output_path;

	std::vector<string> KEYS;
	std::vector<string> SysList;
	std::vector<string> MolList;
	std::vector<string> MonList;
	std::vector<string> LatList;
	std::vector<string> NewtonList;
	std::vector<string> OutputList;
	std::vector<std::string> elems;
	std::vector<string> StateList;
	std::vector<string> ReactionList;


	void PrintList(const std::vector<std::string>&) const;
	std::vector<std::string>& split(const std::string&, char, std::vector<std::string>&) const;
	bool TestNum(std::vector<std::string>&, const string&, int, int, int) const;
	int GetNumStarts(void) const;
	bool CheckParameters(const string&, const string&, int, const std::vector<std::string>&, ParameterStore&) const;
	bool LoadItems(const string&, std::vector<std::string>&, std::vector<std::string>&, std::vector<std::string>&) const;
	bool CheckInput(void);
	bool InSet(const std::vector<std::string>&, const string&) const;
	bool InSet(const vector<int>&, int) const;
	bool InSet(const vector<int>&, int&, int) const;
	bool InSet(const std::vector<std::string>&, int&, const string&) const;
	bool ReadFile(const string&, string&) const;
	bool EvenBrackets(const string&, vector<int>&, vector<int>&) const;
	bool EvenSquareBrackets(const string&, vector<int>&, vector<int>&) const;
	bool MakeLists(int);
	const std::string& GetOutputPath() const;
	std::string ResolvePath(const std::string&) const;

private:
	void parseOutputInfo();
	bool OutputPathExists() const;
};

#endif
