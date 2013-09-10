#include "doc_list_index.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace sdsl;

using idx_type = IDX_TYPE;

const size_t buf_size=1024*128;
char   buffer[buf_size];

template<uint8_t t_width>
struct myline {
    static string parse(char* str) {
        return string(str);
    }
};

template<>
struct myline<0> {
    static vector<uint64_t> parse(char* str) {
        vector<uint64_t> res;
        stringstream ss(str);
        uint64_t x;
        while (ss >> x) {
            res.push_back(x);
        }
        return res;
    }
};

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " index_file pattern_file" << endl;
        cout << " Process all queries with the index." << endl;
        return 1;
    }
    string index_file = string(argv[1]);
    string pattern_file = string(argv[2]);
    idx_type idx;

    using timer = std::chrono::high_resolution_clock;

    std::cout<<"# index_file = "<<index_file<<endl;
    if (!load_from_file(idx, index_file)) {
        std::cerr << "Could not load index file" << std::endl;
        return 1;
    }
    std::cout<<"# pattern_file = "<<pattern_file<<endl;
    std::cout<<"# doc_cnt = "<<idx.doc_cnt()<<endl;
    std::cout<<"# word_cnt = "<<idx.word_cnt()<<endl;
    ifstream in(pattern_file);
    if (!in) {
        std::cerr << "Could not load pattern file" << std::endl;
        return 1;
    }

    using timer = std::chrono::high_resolution_clock;
    size_t q_len = 0;
    size_t q_cnt = 0;
    size_t sum = 0;
    size_t sum_fdt = 0;
    bool tle = false; // flag: time limit exceeded
    auto start = timer::now();
    while (!tle and in.getline(buffer, buf_size)) {
        auto q_start = timer::now();
        typename idx_type::result res;
        auto query = myline<idx_type::WIDTH>::parse(buffer);
        q_len += query.size();
        ++q_cnt;
        size_t x = idx.search(query.begin(), query.end(), res, 10);
        sum += x;
for (auto& r : res) {
            sum_fdt += r.second;
        }
        auto q_time = timer::now()-q_start;
        // single query should not take more then 5 seconds
        if (std::chrono::duration_cast<std::chrono::seconds>(q_time).count() > 5) {
            tle = true;
        }
    }
    auto stop = timer::now();
    auto elapsed = stop-start;
    std::cout<<"# TLE = " << tle << endl;
    std::cout<<"# query_len = "<<q_len/q_cnt<<endl;
    std::cout<<"# queries = " <<q_cnt <<endl;
    std::cout<<"# time_per_query = "<<std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count()/q_cnt <<endl;
    std::cout<<"# check_sum = "<<sum<<endl;
    std::cout<<"# check_sum_fdt = "<<sum_fdt<<endl;
}
