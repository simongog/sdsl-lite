#ifndef SDSL_CONFIG
#define SDSL_CONFIG

#define CMAKE_SOURCE_DIR @CMAKE_SOURCE_DIR@
#cmakedefine divsufsort_FOUND
#cmakedefine divsufsort64_FOUND

#include "uintx_t.hpp"

namespace sdsl{
	namespace constants{ // namespace for library constant
		// size of the buffer for reading and writing data in elements (not in bytes)
		const uint64_t SDSL_BLOCK_SIZE = (uint64_t)1<<22; 
	}
}

#endif
