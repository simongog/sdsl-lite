
#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;

typedef CST_TYPE cst_type;
typedef cst_sada<cst_type::csa_type, lcp_dac<> > fast_cst_type;

int main(int argc, char **argv) {
    if(argc < 3) {
        std::cout << "Usage: " << argv[0] << " input_file tmp_dir" << std::endl;
        std::cout << " Generates a CST and outputs statistics about the distribution of depth." << std::endl;
        std::cout << " Temporary files are stored in tmp_dir." << std::endl;
        return 1;
    }

    const char* input_file = argv[1];
    const char* tmp_dir = argv[2];

    fast_cst_type cst;
    cache_config config(true, tmp_dir);
    construct(cst, input_file, config, fast_cst_type::alphabet_category::WIDTH == 8 ? 1 : 0);

    std::vector<size_t> depth_cnt(100 * 1024 * 1024 + 1, 0);

    for(auto it = cst.begin(); it != cst.end(); ++it) {
        if(it.visit() == 1 && !cst.is_leaf(*it)) {
            const size_t depth = cst.depth(*it);
            depth_cnt[depth]++;
        }
    }

    int last_i = -1;
    size_t sum = 0;
    for(int i = 0; i < 1000000; i++) {
        sum += depth_cnt[i];
        if((i - last_i) > i / 100) {
            std::cout << i << " " << sum << std::endl;
            last_i = i;
        }
    }

    return 0;
}
