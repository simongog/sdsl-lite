#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <random>

using namespace std;
using namespace sdsl;

#ifndef INT_ALPHABET
using csa_t = csa_wt<wt_huff<rrr_vector<63>>>;
uint8_t num_bytes = 1;
#else
using csa_t = csa_wt<wt_int<rrr_vector<63>>>;
uint8_t num_bytes = 0;
#endif

int main(int argc, char* argv[])
{
    if (argc < 7) {
        cout << "Usage: " << argv[0] << " collection_file collection_csa tmp_dir pattern_length pattern_number pattern_file" << endl;
        cout << " Generates an index and stores result in index_file" << endl;
        cout << " Temporary files are stored in tmp_dir." << endl;
        return 1;
    }
    string collection_file = argv[1];
    string id              = util::basename(collection_file);
    string collection_csa  = argv[2];
    string tmp_dir         = argv[3];
    uint64_t pat_len       = stoull(argv[4]);
    uint64_t pat_num       = stoull(argv[5]);
    string pattern_file    = argv[6];

    csa_t csa;
    cache_config cconfig(false, tmp_dir, id);
    if (!load_from_file(csa, collection_csa)) {
        if (num_bytes == 0) {
            int_vector<> v;
            load_from_file(v, collection_file);
            std::cout<<"v.size()="<<v.size()<<std::endl;
        }
        construct(csa, collection_file, cconfig, num_bytes);
        std::cout<<util::demangle2(typeid(csa_t).name())<<std::endl;
        store_to_file(csa, collection_csa);
    }

    // if pat_len < size of CSA - separators - sentinel
    if (pat_len+1 > csa.size() - csa.bwt.rank(csa.size(), 1)) {
        std::cerr<<"pat_len > " << " length of the documents" << std::endl;
        return 1;
    }

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, csa.size()-pat_len);
    auto dice = bind(distribution, rng);

    ofstream out(pattern_file);
    if (!out) {
        std::cerr<<"Could not open file "<<pattern_file<<" for output"<<std::endl;
        return 1;
    }

    uint64_t pat_cnt=0;
    while (pat_cnt < pat_num) {
        uint64_t pos = dice();
        auto pat = extract(csa, pos, pos+pat_len-1);
        bool valid = true;
        for (uint64_t i=0; valid and i < pat.size(); ++i) {
            // if pattern includes separator or newline in byte sequence
            if (pat[i] == 1 or (num_bytes == 1 and pat[i]=='\n')) {
                valid = false;
            }
        }
        if (valid) {
            if (csa_t::alphabet_category::WIDTH == 0 or
                count(csa, pat.begin(), pat.end()) >= 5) {
                out << pat << "\n";
                ++pat_cnt;
            }
        }
    }
}
