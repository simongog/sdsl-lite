#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <string>
#include <cstdio>

using namespace sdsl;
using namespace std;

uint64_t my_rand64()
{
    int_vector<> v(1);
    util::set_random_bits(v);
    return v[0];
}

int main(int argc, char** argv)
{
    if (argc < 7) {
        printf("Usage ./%s input_file tmp_dir pat_len min_cum_occ pat_file occ_file\n",argv[0]);
        printf(" (1) Build an index for the text T stored in input_file.\n");
        printf("     Temporary files are stored in tmp_dir.\n");
        printf(" (2) Set cum_occ=0\n");
        printf(" (3) Select a random position in T and extracts a pattern P\n");
        printf("     of length pat_len from this position.\n");
        printf(" (3) Get x, the number of occurrences of P in T.\n");
        printf(" (4) cum_occ+=x; if cum_occ>=min_cum_occ exit program.\n");
        printf("     Otherwise go to step (3).");
        printf(" (5) Write pattern in Pizza&Chili format into pat_file.");
        printf(" (6) Write #occurrences in into occ_file.");
        return 1;
    }
    uint64_t min_cum_occ = atoll(argv[4]);
    uint64_t pat_len = atoll(argv[3]);
    assert(pat_len > 0);
    vector<string> pattern;
    string file = util::basename(argv[1]);
    ofstream out(argv[5]);
    ofstream occ_out(argv[6]);
//(1)
    csa_wt<> csa;
    cache_config cconfig(false, argv[2], file);
    construct(csa, argv[1], cconfig, 1);
    if (pat_len>=csa.size()) {
        std::cout << "pat_len=" << pat_len << " too long." << std::endl;
        return 1;
    }
// (2)
    uint64_t cum_sum = 0;
    uint64_t number = 0;
    do {
// (3)
        uint64_t start_idx = my_rand64()%(csa.size()-pat_len);
        auto pat = extract(csa, start_idx, start_idx+pat_len-1);
        pattern.push_back(pat);
// (4)
        uint64_t x = count(csa, pat.begin(), pat.end());
        occ_out<<x<<std::endl;
        cum_sum += x;
        ++number;
    } while (cum_sum < min_cum_occ);
// (5)
    out<<"# number="<<number<<" length="<<pat_len<<" file="<<argv[1];
    out<<" forbidden="<<std::endl;
    for (size_t i=0; i<pattern.size(); ++i) {
        out<<pattern[i];
    }
}
