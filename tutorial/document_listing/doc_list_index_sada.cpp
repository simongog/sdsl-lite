#include "doc_list_index_sada.hpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace sdsl;

using idx_type = doc_list_index_sada<>;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " collection_file" << endl;
        return 1;
    }
    string collection_file = string(argv[1]);
    idx_type idx;
    string idx_file = collection_file+".idx";


    using timer = std::chrono::high_resolution_clock;

    if (!load_from_file(idx, idx_file)) {
        cout << "Generate index for " << collection_file << endl;
        {
            auto start = timer::now();
            idx_type temp_idx(collection_file, 1);
            auto stop = timer::now();
            auto elapsed = stop-start;
            std::cout << "construction time = " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << std::endl;
            store_to_file(temp_idx, idx_file);
            ofstream out(idx_file+".html");
            write_structure<HTML_FORMAT>(temp_idx, out);
        }
        load_from_file(idx, idx_file);
    } else {
        cout << "Loaded index from " << collection_file << endl;
    }

    using timer = std::chrono::high_resolution_clock;
    char   buffer[64];
    size_t q_len = 0;
    size_t q_cnt = 0;
    size_t sum = 0;
    auto start = timer::now();
    while (cin.getline(buffer, 64)) {
        typename idx_type::result res;
        string query(buffer);
        if (q_len !=  query.size()) {
            if (q_len == 0) {
                start = timer::now();
            } else {
                auto stop = timer::now();
                auto elapsed = stop-start;
                std::cout<<q_len<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/q_cnt << std::endl;
                start = timer::now();
                q_cnt = 0;
            }
            q_len = query.size();
        }
        ++q_cnt;
        sum += idx.search(query.begin(), query.end(), res, 10);
    }
    auto stop = timer::now();
    auto elapsed = stop-start;
    std::cout<<q_len<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()/q_cnt << std::endl;

    cerr << "sum = " << sum << endl;
}
