#include "sdsl/construct_sa.hpp"

namespace sdsl
{

void construct_sa_se(cache_config& config)
{
    int_vector<8> text;
    load_from_file(text, cache_file_name(conf::KEY_TEXT, config));

    if (text.size() <= 2) {
        // If text is c$ or $ write suffix array [1, 0] or [0]
        int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, config), std::ios::out, 8, 2);
        if (text.size() == 2) {
            sa.push_back(1);
        }
        sa.push_back(0);
    } else {
        _construct_sa_se<int_vector<8>>(text, cache_file_name(conf::KEY_SA, config), 256, 0);
    }
    register_cache_file(conf::KEY_SA, config);
}

}
