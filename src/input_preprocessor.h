#ifndef INPUT_PREPROCESSORxH
#define INPUT_PREPROCESSORxH

#include <nlohmann/json.hpp>
#include <string>

namespace io::input {

// Reads a JSON or legacy input into a normalized document. Legacy inputs are also written as <stem>.input.json.
bool PrepareInput(const std::string&, nlohmann::ordered_json&);

} // namespace io::input

#endif
