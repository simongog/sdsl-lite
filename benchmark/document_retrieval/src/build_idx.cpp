#include "doc_list_index.hpp"
#include <iostream>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace sdsl;

using idx_type = IDX_TYPE;

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " collection_file tmp_dir index_file" << endl;
        cout << " Generates an index and stores result in index_file" << endl;
        cout << " Temporary files are stored in tmp_dir." << endl;
        return 1;
    }
    string collection_file = argv[1];
    string id              = util::basename(collection_file);
    string tmp_dir         = argv[2];
    string idx_file        = argv[3];

    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();

    idx_type idx;
    cache_config cconfig(false, tmp_dir, id);
    construct(idx, collection_file, cconfig, idx_type::WIDTH==8 ? 1 : 0);
    auto stop = timer::now();
    auto elapsed = stop-start;
    std::cout << "construction time = " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << std::endl;
    store_to_file(idx, idx_file);
}
