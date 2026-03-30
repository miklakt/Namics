#ifndef DEBUG_LOG_H
#define DEBUG_LOG_H

#include <iostream>
#include <ostream>
#include <sstream>

extern bool debug;

// Expands to a compiler-specific function signature used in debug output.
#if defined(_MSC_VER)
#define NAMICS_FUNCTION_SIGNATURE __FUNCSIG__
#elif defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER) || defined(__clang__) || defined(__GNUC__)
#define NAMICS_FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#else
#define NAMICS_FUNCTION_SIGNATURE __func__
#endif

/**
 * Emits one debug log line to stderr.
 *
 * Output format:
 *   [DBG] <file>:<line> | <function_signature> | <message>
 *
 * Example emission:
 *   [DBG] LGrad1.cpp:135 | virtual void LGrad1::propagate(Real*, Real*, int, int, int) |  propagate in LGrad1
 *
 * Warning:
 *   This logger is intentionally mutex-free for now.
 *   In multithreaded runs, log lines may interleave; reintroduce synchronization.
 */
template <typename FillFn>
inline void DBG_Emit(const char* file,
                     int line,
                     const char* function_signature,
                     FillFn&& fill) {
	if (!::debug) return;
	const char* base = (file == nullptr) ? "" : file;
	for (const char* p = base; *p != '\0'; ++p) {
		if (*p == '/' || *p == '\\') base = p + 1;
	}
	std::ostringstream payload;
	fill(payload);

	std::cerr << "[DBG] "
	          << base << ':' << line
	          << " | " << (function_signature ? function_signature : "")
	          << " | " << payload.str();
	const std::string s = payload.str();
	if (s.empty() || s.back() != '\n') {
		std::cerr << '\n';
	}
}

/**
 * Convenience macro for stream-style debug messages.
 * Use it as: NAMICS_DBG("text " << value << std::endl);
 */
#define NAMICS_DBG(message_expr)                                                        \
	do {                                                                                \
		DBG_Emit(__FILE__, __LINE__, NAMICS_FUNCTION_SIGNATURE,                         \
		                         [&](std::ostream& _namics_debug_os) {                  \
			                         _namics_debug_os << message_expr;                  \
		                         });                                                     \
	} while (0)

#endif // DEBUG_LOG_H
