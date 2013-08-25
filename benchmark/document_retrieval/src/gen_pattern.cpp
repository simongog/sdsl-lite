#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <random>

using namespace std;
using namespace sdsl;

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

    csa_wt<wt_huff<rrr_vector<63>>> csa;
    cache_config cconfig(false, tmp_dir, id);
    if (!load_from_file(csa, collection_csa)) {
        construct(csa, collection_csa, cconfig, 1);
        store_to_file(csa, collection_csa);
    }

    // if pat_len < size of CSA - separators - sentinel
    if (pat_len+1 > csa.size() - csa.rank_bwt(csa.size(), 1)) {
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
        string pat = extract<string>(csa, pos, pos+pat_len-1);
        if (pat.find_first_of("\n") == string::npos and
            count(csa, pat.begin(), pat.end()) >= 5) {
            out << pat << "\n";
            ++pat_cnt;
        }
    }
}
