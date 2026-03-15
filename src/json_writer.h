#ifndef JSON_WRITER_H
#define JSON_WRITER_H

#include "namics.h"

#include <nlohmann/json.hpp>
#include <string>

namespace io {
namespace json {

class JsonWriter {
public:
	bool WriteProblem(const std::string& filename,
	                  const nlohmann::ordered_json& problem_object,
	                  const nlohmann::ordered_json& metadata_object,
	                  bool append_existing_file,
	                  bool first_problem_of_run) const;
};

} // namespace json
} // namespace io

#endif
