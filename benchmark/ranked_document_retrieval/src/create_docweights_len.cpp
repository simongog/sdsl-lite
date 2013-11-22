#include <iostream>
#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " collection_file docweight_file" << endl;
        cout << " Generates the docweight file from the collection file" << endl;
        return 1;
    }
    string collection_file = argv[1];
    string docweight_file  = argv[2];

    int_vector_buffer<8> in_buf(collection_file,std::ios::in,1024*1024,8,true);
    if (! in_buf.is_open()) {
        cerr << "Can not open collection file.";
        return 1;
    }
    std::ofstream w_fs(docweight_file);
    if (! w_fs.is_open()) {
        cerr << "Can not open docweight file.";
        return 1;
    }

    size_t cur_doc_size = 0;
    for (size_t i=0; i<in_buf.size(); i++) {
        auto sym = in_buf[i];
        if (sym == 1) {
            w_fs << cur_doc_size << "\n";
            cur_doc_size = 0;
        }
        cur_doc_size++;
    }
}
