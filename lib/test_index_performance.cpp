#include "sdsl/test_index_performance.hpp"

namespace sdsl
{

int_vector<64> get_rnd_positions(uint8_t log_s, uint64_t& mask, uint64_t mod, uint64_t seed)
{
    mask = (1<<log_s)-1;
    int_vector<64> rands(1<<log_s ,0);
    util::set_random_bits(rands, seed);
    if (mod>0) {
        util::mod(rands, mod);
    }
    return rands;
}

}// end namespace
