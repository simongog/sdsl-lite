#ifndef SDSL_CONFIG
#define SDSL_CONFIG

#include "uintx_t.hpp"
#include <map>
#include <string>

namespace sdsl
{
namespace conf  // namespace for library constant
{
// size of the buffer for reading and writing data in elements (not in bytes)
const uint64_t SDSL_BLOCK_SIZE = (uint64_t)1<<22;

const char KEY_BWT[] 		= "bwt";
const char KEY_BWT_INT[]	= "bwt_int";
const char KEY_SA[] 		= "sa";
const char KEY_CSA[] 		= "csa";
const char KEY_CST[] 		= "cst";
const char KEY_ISA[] 		= "isa";
const char KEY_TEXT[] 		= "text";
const char KEY_TEXT_INT[] 	= "text_int";
const char KEY_PSI[] 		= "psi";
const char KEY_LCP[] 		= "lcp";
const char KEY_SAMPLE_CHAR[]= "sample_char";
}
typedef uint64_t int_vector_size_type;

typedef std::map<std::string, std::string> tMSS;

enum format_type {JSON_FORMAT, R_FORMAT, HTML_FORMAT};

enum byte_sa_algo_type {LIBDIVSUFSORT, SE_SAIS};

//! Helper class for construction process
struct cache_config {
    bool 		delete_files;   // Flag which indicates if all files which were created
    // during construction should be deleted.
    std::string dir;    		// Directory for temporary files.
    std::string id;     		// Identifier is part of temporary file names. If
    // id is the empty string, then it will be replace
    // a concatenation of PID and a unique ID inside the
    // current process.
    tMSS 		file_map;		// Files stored during the construction process.
    cache_config(bool f_delete_files=true, std::string f_dir="./", std::string f_id="", tMSS f_file_map=tMSS());
};

//! Helper classes to transform width=0 and width=8 to corresponding text key
template<uint8_t width>
struct key_text_trait {
    static const char* KEY_TEXT;
};

//! Helper classes to transform width=0 and width=8 to corresponding bwt key
template<uint8_t width>
struct key_bwt_trait {
    static const char* KEY_BWT;
};
}

#endif
