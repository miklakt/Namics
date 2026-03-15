#ifndef INPUTxH
#define INPUTxH
#include "namics.h"

template <typename T>
inline bool ParseStrict(const std::string& s, T& value) {
	if (s.empty()) return false;
	std::istringstream stream(s);
	stream >> std::noskipws >> value;
	return !stream.fail() && stream.eof();
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

class Input {
public:
	Input(const std::string&);

	~Input();

	std::string json_path;
	bool Input_error;
	std::string output_path;
	ParameterStore starts;
	int active_start;

	std::vector<std::string> SysList;
	std::vector<std::string> MolList;
	std::vector<std::string> MonList;
	std::vector<std::string> LatList;
	std::vector<std::string> NewtonList;
	std::vector<std::string> OutputList;
	std::vector<std::string> StateList;
	std::vector<std::string> ReactionList;


	std::vector<std::string>& split(const std::string&, char, std::vector<std::string>&) const;
	int GetNumStarts(void) const;
	const ParameterStore& Parameters(const std::string&, const std::string&, int) const;
	ParameterStore LoadItems(const std::string&) const;
	bool CheckInput(void);
	bool EvenBrackets(const std::string&, std::vector<int>&, std::vector<int>&) const;
	bool EvenSquareBrackets(const std::string&, std::vector<int>&, std::vector<int>&) const;
	bool MakeLists(int);
	const std::string& GetOutputPath() const;
	std::string ResolvePath(const std::string&) const;
	const ParameterStore& operator[](const std::string&) const;
	const ParameterStore& Start(int) const;

private:
	void UpdateOutputPath();
	bool OutputPathExists() const;
};

#endif
