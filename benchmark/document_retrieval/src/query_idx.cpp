#include "doc_list_index.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace sdsl;

using idx_type = IDX_TYPE;

const size_t buf_size=1024;
char   buffer[buf_size];

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
    auto start = timer::now();
    while (in.getline(buffer, buf_size)) {
        typename idx_type::result res;
        string query(buffer);
        q_len += query.size();
        ++q_cnt;
        size_t x = idx.search(query.begin(), query.end(), res, 10);
        sum += x;
        for (auto& r : res) {
            sum_fdt += r.second;
        }
    }
    auto stop = timer::now();
    auto elapsed = stop-start;
    std::cout<<"# query_len = "<<q_len/q_cnt<<std::endl;
    std::cout<<"# queries = " <<q_cnt << std::endl;
    std::cout<<"# time_per_query = "<<std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count()/q_cnt << std::endl;
    std::cout<<"# check_sum = "<<sum<<std::endl;
    std::cout<<"# check_sum_fdt = "<<sum_fdt<<std::endl;
}
