#ifndef IO_UTILSxH
#define IO_UTILSxH

#include "namics.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace io {

inline std::ofstream OpenOutputStream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out) {
	return std::ofstream(filename, mode);
}

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

inline bool ReadSanitizedLines(const std::string& filename, std::vector<std::string>& lines) {
	std::string content;
	if (!ReadSanitizedFile(filename, content)) return false;
	lines.clear();
	std::stringstream ss(content);
	std::string item;
	while (std::getline(ss, item, '#')) {
		if (!item.empty()) lines.push_back(item);
	}
	return true;
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
	std::ifstream input(filename.c_str());
	if (!input.is_open()) {
		std::cout << "inputfile " << filename << " is not found. Read guess for initial guess failed" << std::endl;
		return false;
	}

	std::string charge_flag;
	int mon_count = 0;
	int state_count = 0;
	if (!(input >> method >> mx >> my >> mz >> fjc >> charge_flag >> mon_count)) {
		return false;
	}
	charged = (charge_flag == "true");

	monlist.clear();
	monlist.reserve(std::max(0, mon_count));
	for (int i = 0; i < mon_count; ++i) {
		std::string name;
		if (!(input >> name)) return false;
		monlist.push_back(name);
	}

	if (!(input >> state_count)) return false;
	statelist.clear();
	statelist.reserve(std::max(0, state_count));
	for (int i = 0; i < state_count; ++i) {
		std::string name;
		if (!(input >> name)) return false;
		statelist.push_back(name);
	}

	if (readx == 0) return true;

	int m = 0;
	if (my == 0) {
		m = (mx + 2 * fjc);
	} else if (mz == 0) {
		m = (mx + 2 * fjc) * (my + 2 * fjc);
	} else {
		m = (mx + 2 * fjc) * (my + 2 * fjc) * (mz + 2 * fjc);
	}
	int total = (mon_count + state_count) * m;
	if (charged) total += m;
	for (int i = 0; i < total; ++i) {
		if (!(input >> x[i])) return false;
	}

	return true;
}

template <typename RealT>
inline bool WriteInitialGuess(const std::string& filename,
                              const std::string& method,
                              int mx,
                              int my,
                              int mz,
                              int fjc,
                              bool charged,
                              const std::vector<std::string>& monlist,
                              const std::vector<std::string>& statelist,
                              const RealT* values,
                              int value_count) {
	std::ofstream out = OpenOutputStream(filename, std::ios::out | std::ios::trunc);
	if (!out.is_open()) return false;

	out << method << '\n';
	out << ' ' << mx << '\t' << my << '\t' << mz << '\t' << fjc << '\n';
	out << (charged ? "true" : "false") << '\n';
	out << static_cast<int>(monlist.size()) << '\n';
	for (const auto& mon : monlist) out << mon << '\n';
	out << static_cast<int>(statelist.size()) << '\n';
	for (const auto& state : statelist) out << state << '\n';

	out << std::setprecision(std::numeric_limits<RealT>::digits10 + 2);
	for (int i = 0; i < value_count; ++i) out << values[i] << '\n';

	return static_cast<bool>(out);
}

} // namespace io

#endif
