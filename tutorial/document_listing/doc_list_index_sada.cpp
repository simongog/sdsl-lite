#include "doc_list_index_sada.hpp"
#include <iostream>

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

    if (!load_from_file(idx, idx_file)) {
        cout << "Generate index for " << collection_file << endl;
        {
            idx_type temp_idx(collection_file, 1);
            store_to_file(temp_idx, idx_file);
        }
        load_from_file(idx, idx_file);
    } else {
        cout << "Loaded index from " << collection_file << endl;
    }

    string query;
    while (cin >> query) {
        typename idx_type::result res;
        uint64_t docs = idx.search(query.begin(), query.end(), res);
        cout << "Found " << docs << " documents." << endl;
    }
}
