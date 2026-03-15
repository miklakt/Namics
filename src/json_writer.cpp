#include "json_writer.h"

#include <fstream>
#include <iomanip>

namespace io {
namespace json {

bool JsonWriter::WriteProblem(const std::string& filename,
                              const nlohmann::ordered_json& problem_object,
                              const nlohmann::ordered_json& metadata_object,
                              bool append_existing_file,
                              bool first_problem_of_run) const {
	nlohmann::ordered_json document = {
		{"problems", nlohmann::ordered_json::array({problem_object})},
		{"metadata", metadata_object}
	};
	if (!(first_problem_of_run && !append_existing_file)) {
		std::ifstream in(filename.c_str());
		if (in.is_open()) {
			try {
				nlohmann::ordered_json existing;
				in >> existing;
				if (!first_problem_of_run) {
					if (existing.is_object() && existing.contains("problems") && existing["problems"].is_array()) {
						existing["problems"].push_back(problem_object);
						document = std::move(existing);
					} else if (existing.is_array() && !existing.empty() && existing.back().is_object() &&
					           existing.back().contains("problems") && existing.back()["problems"].is_array()) {
						existing.back()["problems"].push_back(problem_object);
						document = std::move(existing);
					}
				} else if (append_existing_file) {
					if (existing.is_array()) {
						existing.push_back(document);
						document = std::move(existing);
					} else if (existing.is_object() && existing.contains("problems") && existing["problems"].is_array()) {
						document = nlohmann::ordered_json::array({existing, document});
					}
				}
			} catch (...) {
			}
		}
	}

	std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc);
	if (!out.is_open()) return false;
	out << std::setw(2) << document << '\n';
	return static_cast<bool>(out);
}

} // namespace json
} // namespace io
