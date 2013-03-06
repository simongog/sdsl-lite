#ifndef INCLUDED_SDSL_LCP_CONSTRUCT_HELPER
#define INCLUDED_SDSL_LCP_CONSTRUCT_HELPER

#include "sdsl/int_vector.hpp"

namespace sdsl{


void insert_lcp_values(int_vector<> &partial_lcp, bit_vector &index_done, std::string lcp_file, uint64_t max_lcp_value, uint64_t lcp_value_offset);

template<class tWT>
void create_C_array(std::vector<uint64_t> &C, const tWT &wt){ 
	uint64_t quantity;                          // quantity of characters in interval
	std::vector<unsigned char> cs(wt.sigma);      // list of characters in the interval
	std::vector<uint64_t> rank_c_i(wt.sigma);    // number of occurrence of character in [0 .. i-1]
	std::vector<uint64_t> rank_c_j(wt.sigma);    // number of occurrence of character in [0 .. j-1]

    C = std::vector<uint64_t>(257, 0);				
    wt.interval_symbols(0, wt.size(), quantity, cs, rank_c_i, rank_c_j);
    for(uint64_t i=0; i<quantity; ++i) {
        unsigned char c = cs[i];
        C[c+1] = rank_c_j[i];
    }
    for(uint64_t i=1; i<C.size()-1; ++i) {
        C[i+1] += C[i];
    }
}

}

#endif
