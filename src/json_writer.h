#ifndef JSON_WRITER_H
#define JSON_WRITER_H

#include "namics.h"

#include <memory>
#include <string>

namespace io {
namespace json {

class JsonWriter {
public:
	bool WriteProblem(const std::string& filename,
	                  const std::string& problem_object,
	                  const std::string& metadata_object,
	                  bool append_existing_file,
	                  bool first_problem_of_run) const;
};

inline std::shared_ptr<JsonWriter> SharedJsonWriter() {
	static std::shared_ptr<JsonWriter> instance = std::make_shared<JsonWriter>();
	return instance;
}

} // namespace json
} // namespace io

#endif
