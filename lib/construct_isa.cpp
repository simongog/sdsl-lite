/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog
*/
#include "sdsl/construct_isa.hpp"
#include <string>

namespace sdsl
{

void construct_isa(cache_config& config)
{
    typedef int_vector<>::size_type size_type;
    if (!cache_file_exists(conf::KEY_ISA, config)) {   // if isa is not already on disk => calculate it
        int_vector_buffer<> sa_buf(cache_file_name(conf::KEY_SA, config));
        if (!sa_buf.is_open()) {
            throw std::ios_base::failure("cst_construct: Cannot load SA from file system!");
        }
        int_vector<> isa(sa_buf.size());
        for (size_type i=0; i < isa.size(); ++i) {
            isa[ sa_buf[i] ] = i;
        }
        store_to_cache(isa, conf::KEY_ISA, config);
    }
}

}// end namespace
