#ifndef DEBUG_LOG_H
#define DEBUG_LOG_H

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string>
#include <typeinfo>
#if defined(__GNUG__)
#include <cxxabi.h>
#endif

extern bool debug;

#if defined(_MSC_VER)
#define NAMICS_FUNCTION_SIGNATURE __FUNCSIG__
#elif defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) || defined(__clang__) || defined(__GNUC__)
#define NAMICS_FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#else
#define NAMICS_FUNCTION_SIGNATURE __func__
#endif

namespace namics::debuglog {

inline const char* Basename(const char* path) {
	if (path == nullptr) return "";
	const char* file = path;
	for (const char* p = path; *p != '\0'; ++p) {
		if (*p == '/' || *p == '\\') file = p + 1;
	}
	return file;
}

inline std::mutex& OutputMutex() {
	static std::mutex mutex;
	return mutex;
}

inline std::string Demangle(const char* name) {
	if (name == nullptr) return "";
#if defined(__GNUG__)
	int status = 0;
	char* demangled = abi::__cxa_demangle(name, nullptr, nullptr, &status);
	if (status != 0 || demangled == nullptr) return std::string(name);
	std::string result(demangled);
	std::free(demangled);
	return result;
#else
	return std::string(name);
#endif
}

template <typename FillFn>
inline void Emit(const char* file,
                 int line,
                 const char* function_signature,
                 const void* object_ptr,
                 const char* object_type_name,
                 FillFn&& fill) {
	if (!debug) return;
	std::ostringstream payload;
	fill(payload);

	std::lock_guard<std::mutex> lock(OutputMutex());
	std::cerr << "[DBG] "
	          << Basename(file) << ':' << line
	          << " | " << (function_signature ? function_signature : "");
	if (object_ptr != nullptr) {
		std::cerr << " | this=" << object_ptr;
		if (object_type_name != nullptr && object_type_name[0] != '\0') {
			std::cerr << " type=" << Demangle(object_type_name);
		}
	}
	std::cerr << " | " << payload.str();
	const std::string s = payload.str();
	if (s.empty() || s.back() != '\n') {
		std::cerr << '\n';
	}
}

} // namespace namics::debuglog

#define NAMICS_DBG(message_expr)                                                           \
	do {                                                                                   \
		::namics::debuglog::Emit(__FILE__, __LINE__, NAMICS_FUNCTION_SIGNATURE, nullptr,  \
		                         nullptr,                                                   \
		                         [&](std::ostream& _namics_debug_os) {                     \
			                         _namics_debug_os << message_expr;                     \
		                         });                                                        \
	} while (0)

#define NAMICS_DBG_THIS(message_expr)                                                      \
	do {                                                                                   \
		::namics::debuglog::Emit(__FILE__, __LINE__, NAMICS_FUNCTION_SIGNATURE, this,      \
		                         typeid(*this).name(),                                      \
		                         [&](std::ostream& _namics_debug_os) {                     \
			                         _namics_debug_os << message_expr;                     \
		                         });                                                        \
	} while (0)

#endif // DEBUG_LOG_H
