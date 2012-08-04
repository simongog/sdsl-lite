#include <sdsl/wt_huff.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "usage: "<<argv[0]<< " file_name" << std::endl;
    } else {
        int_vector_file_buffer<8> file_buf;
        file_buf.load_from_plain(argv[1]);
        wt_huff<> wt(file_buf, file_buf.int_vector_size);

        for (size_t i=0; i<wt.size() ; ++i) {
            cout << wt[i];
        }
        cout << endl;
    }
}
