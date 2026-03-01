#ifndef INPUTxH
#define INPUTxH
#include "namics.h"
#include "output_info.h"
#include <functional>

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
	OutputInfo output_info;

	std::vector<string> KEYS;
	std::vector<string> SysList;
	std::vector<string> MolList;
	std::vector<string> MonList;
	std::vector<string> LatList;
	std::vector<string> AliasList;
	std::vector<string> NewtonList;
	std::vector<string> OutputList;
	std::vector<string> MicroList;
	std::vector<string> VarList;
	std::vector<std::string> elems;
	std::vector<string> StateList;
	std::vector<string> ReactionList;


	void PrintList(const std::vector<std::string>&) const;
	std::vector<std::string>& split(const std::string&, char, std::vector<std::string>&) const;
	bool IsDigit(const string&) const;
	int Get_int(const string&, int) const;
	bool Get_int(const string&, int&, const std::string&) const;
	bool Get_int(const string&, int&, int, int, const std::string&) const;
	string Get_string(const string&, const string&) const;
	bool Get_string(const string&, string&, const std::string&) const;
	bool Get_string(const string&, string&, const std::vector<std::string>&, const std::string&) const;
	Real Get_Real(const string&, Real) const;
	bool Get_Real(const string&, Real&, const std::string&) const;
	bool Get_Real(const string&, Real&, Real, Real, const std::string&) const;
	bool Get_bool(const string&, bool) const;
	bool Get_bool(const string&, bool&, const std::string&) const;
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
	bool ArePair(char, char) const;
	bool EvenBrackets(const string&, vector<int>&, vector<int>&) const;
	bool EvenSquareBrackets(const string&, vector<int>&, vector<int>&) const;
	bool MakeLists(int);

private:
	void parseOutputInfo();
};

#endif
