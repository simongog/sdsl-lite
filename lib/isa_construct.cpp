/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog
*/
#include "sdsl/isa_construct.hpp"
#include <string>

namespace sdsl
{

void construct_isa(cache_config &config) {
    typedef int_vector<>::size_type size_type;
	tMSS::const_iterator key = config.file_map.find(constants::KEY_ISA);
    if ( config.file_map.end() == key ) { // if isa is not already on disk => calculate it
        write_R_output("cst", "construct ISA", "begin", 1, 0);
        int_vector<> isa;
        if (!util::load_from_file(isa, config.file_map[constants::KEY_SA].c_str())) {
            throw std::ios_base::failure("cst_construct: Cannot load SA from file system!");
        }
        {
            int_vector_file_buffer<> sa_buf(config.file_map[constants::KEY_SA].c_str());

            for (size_type i=0, r_sum=0, r = sa_buf.load_next_block(); r_sum < isa.size();) {
                for (; i<r_sum+r; ++i) {
                    isa[ sa_buf[i-r_sum] ] = i;
                }
                r_sum += r; r = sa_buf.load_next_block();
            }
        }
		util::store_to_cache(isa, constants::KEY_ISA, config);
        write_R_output("cst", "construct ISA", "end", 1, 0);
    }
}

}// end namespace
