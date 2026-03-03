#include "json_writer.h"

#include <cctype>
#include <fstream>
#include <sstream>

namespace io {
namespace json {
namespace {

std::string TrimCopy(const std::string& s) {
	size_t begin = 0;
	while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin]))) ++begin;
	size_t end = s.size();
	while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1]))) --end;
	return s.substr(begin, end - begin);
}

bool WriteText(const std::string& filename, const std::string& text) {
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc);
	if (!out.is_open()) return false;
	out << text;
	return static_cast<bool>(out);
}

std::string BuildDocument(const std::string& problems_content, const std::string& metadata_object) {
	std::string doc;
	doc.reserve(problems_content.size() + metadata_object.size() + 64);
	doc.append("{\n");
	doc.append("  \"problems\": [\n");
	doc.append(problems_content);
	doc.append("\n  ],\n");
	doc.append("  \"metadata\": ");
	doc.append(metadata_object);
	doc.append("\n}");
	return doc;
}

size_t FindMatchingBracket(const std::string& text, size_t open_index) {
	if (open_index >= text.size() || text[open_index] != '[') return std::string::npos;
	int depth = 0;
	bool in_string = false;
	bool escape = false;
	for (size_t i = open_index; i < text.size(); ++i) {
		const char c = text[i];
		if (in_string) {
			if (escape) {
				escape = false;
			} else if (c == '\\') {
				escape = true;
			} else if (c == '"') {
				in_string = false;
			}
			continue;
		}
		if (c == '"') {
			in_string = true;
			continue;
		}
		if (c == '[') {
			++depth;
			continue;
		}
		if (c == ']') {
			--depth;
			if (depth == 0) return i;
		}
	}
	return std::string::npos;
}

bool HasNonWhitespace(const std::string& text, size_t from, size_t to_exclusive) {
	for (size_t i = from; i < to_exclusive; ++i) {
		if (!std::isspace(static_cast<unsigned char>(text[i]))) return true;
	}
	return false;
}

bool ContainsLegacyProblemPayload(const std::string& text, size_t from, size_t to_exclusive) {
	if (from >= to_exclusive || to_exclusive > text.size()) return false;
	const std::string payload = text.substr(from, to_exclusive - from);
	const bool has_profiles_rows = payload.find("\"profiles\"") != std::string::npos &&
	                              payload.find("\"rows\"") != std::string::npos;
	const bool has_scalars = payload.find("\"scalars\"") != std::string::npos;
	return has_profiles_rows || has_scalars;
}

bool AppendProblemToDocumentObject(const std::string& document_object,
                                   const std::string& problem_object,
                                   std::string& updated) {
	if (document_object.empty() || document_object.front() != '{' || document_object.back() != '}') return false;
	const size_t problems_key = document_object.find("\"problems\"");
	if (problems_key == std::string::npos) return false;
	const size_t open = document_object.find('[', problems_key);
	if (open == std::string::npos) return false;
	const size_t close = FindMatchingBracket(document_object, open);
	if (close == std::string::npos) return false;
	if (ContainsLegacyProblemPayload(document_object, open + 1, close)) return false;
	const bool has_entries = HasNonWhitespace(document_object, open + 1, close);
	updated = document_object.substr(0, close);
	if (has_entries) updated.append(",\n");
	else updated.push_back('\n');
	updated.append(problem_object);
	updated.append(document_object.substr(close));
	return true;
}

bool AppendObjectToRootArray(const std::string& root_array,
                             const std::string& object,
                             std::string& updated) {
	if (root_array.empty() || root_array.front() != '[') return false;
	const size_t close = FindMatchingBracket(root_array, 0);
	if (close == std::string::npos) return false;
	const bool has_entries = HasNonWhitespace(root_array, 1, close);
	updated = root_array.substr(0, close);
	if (has_entries) updated.append(",\n");
	else updated.push_back('\n');
	updated.append(object);
	updated.append(root_array.substr(close));
	return true;
}

bool AppendProblemToLastDocumentInRootArray(const std::string& root_array,
                                            const std::string& problem_object,
                                            std::string& updated) {
	if (root_array.empty() || root_array.front() != '[' || root_array.back() != ']') return false;
	const size_t problems_key = root_array.rfind("\"problems\"");
	if (problems_key == std::string::npos) return false;
	const size_t open = root_array.find('[', problems_key);
	if (open == std::string::npos) return false;
	const size_t close = FindMatchingBracket(root_array, open);
	if (close == std::string::npos) return false;
	if (ContainsLegacyProblemPayload(root_array, open + 1, close)) return false;
	const bool has_entries = HasNonWhitespace(root_array, open + 1, close);
	updated = root_array.substr(0, close);
	if (has_entries) updated.append(",\n");
	else updated.push_back('\n');
	updated.append(problem_object);
	updated.append(root_array.substr(close));
	return true;
}

} // namespace

bool JsonWriter::WriteProblem(const std::string& filename,
                              const std::string& problem_object,
                              const std::string& metadata_object,
                              bool append_existing_file,
                              bool first_problem_of_run) const {
	const std::string document = BuildDocument(problem_object, metadata_object);
	if (first_problem_of_run && !append_existing_file) {
		return WriteText(filename, document + "\n");
	}

	std::ifstream in(filename.c_str());
	if (!in.is_open()) {
		return WriteText(filename, document + "\n");
	}

	std::stringstream buffer;
	buffer << in.rdbuf();
	const std::string trimmed = TrimCopy(buffer.str());
	if (trimmed.empty()) {
		return WriteText(filename, document + "\n");
	}

	if (!first_problem_of_run) {
		if (trimmed.front() == '{' && trimmed.back() == '}') {
			std::string updated;
			if (AppendProblemToDocumentObject(trimmed, problem_object, updated)) {
				return WriteText(filename, updated + "\n");
			}
			return WriteText(filename, document + "\n");
		}
		if (trimmed.front() == '[' && trimmed.back() == ']') {
			std::string updated;
			if (AppendProblemToLastDocumentInRootArray(trimmed, problem_object, updated)) {
				return WriteText(filename, updated + "\n");
			}
			return WriteText(filename, document + "\n");
		}
		return WriteText(filename, document + "\n");
	}

	if (trimmed.front() == '[' && trimmed.back() == ']') {
		std::string updated;
		if (AppendObjectToRootArray(trimmed, document, updated)) {
			return WriteText(filename, updated + "\n");
		}
		return WriteText(filename, document + "\n");
	}

	if (trimmed.front() == '{' && trimmed.back() == '}') {
		const size_t problems_key = trimmed.find("\"problems\"");
		if (problems_key != std::string::npos) {
			const size_t open = trimmed.find('[', problems_key);
			if (open != std::string::npos) {
				const size_t close = FindMatchingBracket(trimmed, open);
				if (close != std::string::npos && ContainsLegacyProblemPayload(trimmed, open + 1, close)) {
					return WriteText(filename, document + "\n");
				}
			}
		}
		if (problems_key == std::string::npos) {
			return WriteText(filename, document + "\n");
		}
		std::string updated = "[\n";
		updated.append(trimmed);
		updated.append(",\n");
		updated.append(document);
		updated.append("\n]\n");
		return WriteText(filename, updated);
	}

	return WriteText(filename, document + "\n");
}

} // namespace json
} // namespace io
