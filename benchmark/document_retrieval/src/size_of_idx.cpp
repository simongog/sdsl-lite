#include "doc_list_index.hpp"
#include <iostream>
#include <fstream>

using namespace std;
using namespace sdsl;

using idx_type = IDX_TYPE;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " index_file [html_file]" << endl;
        return 1;
    }
    string index_file = string(argv[1]);
    idx_type idx;

    if (!load_from_file(idx, index_file)) {
        std::cerr << "Could not load index file" << std::endl;
        return 1;
    }
    if (argc > 2) {
        ofstream out(argv[2]);
        write_structure<HTML_FORMAT>(idx, out);
    } else {
        cout<<"# index_file = "<<index_file<<endl;
        cout << "# size = " << size_in_bytes(idx) << endl;
    }
}
