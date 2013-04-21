#include "types.hpp"
#include <cstdlib>
#include <iostream>
#include <queue>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "Usage: ./" << argv[0] << " file tmp_dir" << endl;
        return 1;
    }
    cache_config cconfig(false, argv[2], util::basename(argv[1]));
    if (!cache_file_exists(constants::KEY_LCP, cconfig)) {
        // construct cst to get LCP array
        cst_t cst;
        construct(cst, argv[1], cconfig, 1);
    }
    vector<uint64_t> freq;
    string lcp_file = cache_file_name(constants::KEY_LCP, cconfig);
    int_vector_file_buffer<> lcp_buf(lcp_file);
    for (uint64_t i=0, r_sum=0, r=0; r_sum < lcp_buf.int_vector_size;) {
        for (; i < r_sum+r; ++i) {
            uint64_t x = lcp_buf[i-r_sum];
            if (x >= freq.size()) {
                freq.insert(freq.end(), x+1-freq.size(), 0);
            }
            ++freq[x];
        }
        r_sum += r; r=lcp_buf.load_next_block();
    }
    uint64_t check=0;
    for (size_t i=0; i<freq.size(); ++i) {
        check+=freq[i];
        cout << freq[i] << endl;
    }
    return check != lcp_buf.int_vector_size;
}
