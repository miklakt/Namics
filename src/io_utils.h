#ifndef IO_UTILSxH
#define IO_UTILSxH

#include "namics.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace io {

class Writer {
public:
	FILE* OpenRaw(const std::string& filename, const char* mode) const {
		return std::fopen(filename.c_str(), mode);
	}

	void Close(FILE* file) const {
		if (file != nullptr) std::fclose(file);
	}

	void Flush(FILE* file) const {
		if (file != nullptr) std::fflush(file);
	}

	bool Exists(const std::string& filename) const {
		std::ifstream file(filename.c_str());
		return static_cast<bool>(file);
	}

	std::ofstream OpenOutputStream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out) const {
		return std::ofstream(filename, mode);
	}

	template <typename... Args>
	int Writef(FILE* file, const char* format, Args... args) const {
		if (file == nullptr) return -1;
		if constexpr (sizeof...(Args) == 0) {
			return std::fputs(format, file);
		} else {
			return std::fprintf(file, format, args...);
		}
	}

	bool WriteLineEnding(FILE* file, bool dos) const {
		return Writef(file, dos ? "\r\n" : "\n") >= 0;
	}

	template <typename T>
	bool WriteScalar(FILE* file, const char* format, const T& value, bool dos = false) const {
		if (Writef(file, format, value) < 0) return false;
		return WriteLineEnding(file, dos);
	}

	template <typename T>
	bool WriteVector(FILE* file, std::span<const T> values, const char* format, const char* separator = "\t", bool dos = false) const {
		for (size_t i = 0; i < values.size(); ++i) {
			if (Writef(file, format, values[i]) < 0) return false;
			if (i + 1 < values.size() && Writef(file, "%s", separator) < 0) return false;
		}
		return WriteLineEnding(file, dos);
	}

	template <typename RealT>
	bool WriteInitialGuess(const std::string& filename,
	                      const std::string& method,
	                      int mx,
	                      int my,
	                      int mz,
	                      int fjc,
	                      bool charged,
	                      const std::vector<std::string>& monlist,
	                      const std::vector<std::string>& statelist,
	                      const RealT* values,
	                      int value_count) const {
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
};

class RangeReader {
public:
	bool ReadSanitizedFile(const std::string& filename, std::string& buffer) const {
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

	bool ReadSanitizedLines(const std::string& filename, std::vector<std::string>& lines) const {
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
};

class InitialGuessReader {
public:
	template <typename RealT>
	bool ReadInitialGuess(const std::string& filename,
	                     RealT* x,
	                     std::string& method,
	                     std::vector<std::string>& monlist,
	                     std::vector<std::string>& statelist,
	                     bool& charged,
	                     int& mx,
	                     int& my,
	                     int& mz,
	                     int& fjc,
	                     int readx) const {
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
};

inline std::shared_ptr<Writer> SharedWriter() {
	static std::shared_ptr<Writer> instance = std::make_shared<Writer>();
	return instance;
}

inline std::shared_ptr<RangeReader> SharedRangeReader() {
	static std::shared_ptr<RangeReader> instance = std::make_shared<RangeReader>();
	return instance;
}

inline std::shared_ptr<InitialGuessReader> SharedInitialGuessReader() {
	static std::shared_ptr<InitialGuessReader> instance = std::make_shared<InitialGuessReader>();
	return instance;
}

} // namespace io

#endif
