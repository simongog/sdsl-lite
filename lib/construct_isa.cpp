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
    tMSS::const_iterator key = config.file_map.find(conf::KEY_ISA);
    if (config.file_map.end() == key) {   // if isa is not already on disk => calculate it
        int_vector<> isa;
        if (!load_from_file(isa, config.file_map[conf::KEY_SA])) {
            throw std::ios_base::failure("cst_construct: Cannot load SA from file system!");
        }
        {
            int_vector_buffer<> sa_buf(config.file_map[conf::KEY_SA]);

            for (size_type i=0; i < isa.size(); ++i) {
                isa[ sa_buf[i] ] = i;
            }
        }
        store_to_cache(isa, conf::KEY_ISA, config);
    }
}

}// end namespace
