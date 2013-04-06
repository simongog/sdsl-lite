#include "sdsl/int_vector.hpp"
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 5) {
        cout << "Usage: " << argv[0] << " FILE SIZE WIDTH DEFAULT_VALUE" << endl;
        cout << " (1) Generates an int_vector<>(SIZE, DEFAULT_VALUE, WIDTH)" << endl;
        cout << "     Vector will be initialized with random bits, if " << endl;
        cout << "     DEFAULT_VALUE is r." << endl;
        cout << " (2) Stores the vector to FILE." << endl;
        return 1;
    }

    uint64_t size  = atoll(argv[2]);
    uint64_t width = atoll(argv[3]);
    if ('r' == argv[4][1]) {
        uint64_t default_value = atoll(argv[4]);
        return store_to_file(int_vector<>(size, default_value, width), argv[1]);
    } else {
        int_vector<> iv(size, 0, width);
        util::set_random_bits(iv);
        store_to_file(iv, argv[1]);
    }
}
