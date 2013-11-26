#include "doc_rank_index.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace sdsl;
using namespace std::chrono;

using idx_type = IDX_TYPE;

const size_t buf_size = 1024 * 1024;
char buffer[buf_size];

template <uint8_t t_width>
struct myline {
    static vector<string> parse(char* str) {
        std::string line(str);
        std::vector<std::string> query_terms;
        std::string::size_type p = 0;
        std::string::size_type q;
        while ((q = line.find(',', p)) != std::string::npos) {
            query_terms.emplace_back(line, p, q - p);
            p = q + 1;
        }
        if (p != line.size()) {
            query_terms.emplace_back(line, p);
        }
        return query_terms;
    }
};

template <>
struct myline<0> {
    static vector<vector<uint64_t>> parse(char* str) {
        vector<vector<uint64_t>> query_terms;
        vector<uint64_t> res;
        stringstream ss(str);
        uint64_t x;
        while (ss >> x) {
            if (x == 1) {
                query_terms.push_back(res);
                res.clear();
            }
            res.push_back(x);
        }
        return query_terms;
    }
};

int main(int argc, char* argv[])
{
    using timer = std::chrono::high_resolution_clock;
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " index_file pattern_file k"
                  << endl;
        std::cerr << " Process all queries with the index." << endl;
        return 1;
    }
    string index_file = string(argv[1]);
    string pattern_file = string(argv[2]);
    uint64_t k = std::strtoul(argv[3], NULL, 10);

    idx_type idx;
    if (!load_from_file(idx, index_file)) {
        std::cerr << "Could not load index file" << std::endl;
        return 1;
    }
    ifstream in(pattern_file);
    if (!in) {
        std::cerr << "Could not load pattern file" << std::endl;
        return 1;
    }
    size_t qid = 1;
    while (in.getline(buffer, buf_size)) {
        auto query = myline<idx_type::WIDTH>::parse(buffer);
        auto q_start = timer::now();
        auto res = idx.search(query, k);
        auto elapsed = timer::now() - q_start;

        std::cout << qid << ";" << index_file << ";" << k << ";" << pattern_file
                  << ";" << duration_cast<microseconds>(elapsed).count() << ";"
                  << query.size() << ";" << res.total_range_sizes()
                  << std::endl;
        qid++;
    }
}
