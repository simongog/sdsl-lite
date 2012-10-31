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

		const char KEY_BWT[] 	= "bwt";
		const char KEY_SA[] 	= "sa";
		const char KEY_ISA[] 	= "isa";
		const char KEY_TEXT[] 	= "text";
		const char KEY_PSI[] 	= "psi";
		const char KEY_LCP[] 	= "lcp";
	}
	typedef uint64_t int_vector_size_type;
}

#endif
