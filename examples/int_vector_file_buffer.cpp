#include <sdsl/int_vector.hpp>
#include <string>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " temp_file" << endl;
        cout << "(1) Generates an int_vector" << endl;
        cout << "(2) Stores it to temp_file" << endl;
        cout << "(3) Streams the content of temp_file" << endl;
        cout << "(4) Deletes temp_file" << endl;
        return 1;
    }
    string tmp_file = argv[1];
    size_t size = 10000000;
    // create an int_vector
    int_vector<> v(size);

    // initialize vector with random bits, seed of the random process = 42
    util::set_random_bits(v, 42);

    // store int_vector to temporary file
    store_to_file(v, tmp_file);

    // open file buffer
    int_vector<>::file_buffer v_buf(tmp_file);

    // stream vector data from disk and compare it with in-memory data
    for (size_t i = 0, r_sum = 0, r = v_buf.load_next_block(); i < v.size();) {
        for (; i < r_sum +r; ++i) {
            if (v[i] != v_buf[i-r_sum]) {
                std::cerr << "ERROR: v["<< i << "] != v_buf[" << i-r_sum << "]" << std::endl;
                return 1;
            }
        }
        r_sum += r; r = v_buf.load_next_block();
    }
    // remove temporary file
    sdsl::remove(tmp_file);
}
