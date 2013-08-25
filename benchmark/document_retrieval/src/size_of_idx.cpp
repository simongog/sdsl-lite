#include "doc_list_index.hpp"
#include <iostream>

using namespace std;
using namespace sdsl;

using idx_type = IDX_TYPE;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " index_file" << endl;
        return 1;
    }
    string index_file = string(argv[1]);
    idx_type idx;

    cout<<"# index_file = "<<index_file<<endl;
    if (!load_from_file(idx, index_file)) {
        std::cerr << "Could not load index file" << std::endl;
        return 1;
    }
    cout << "# size = " << size_in_bytes(idx) << endl;
}
