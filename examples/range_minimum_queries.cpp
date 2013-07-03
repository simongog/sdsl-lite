#include <sdsl/rmq_support.hpp> // include header for range minimum queries
#include <iostream>
#include <cstdlib>

using namespace sdsl;
using namespace std;

int main(int argc, char** argv)
{
    uint64_t len = 50;
    if (argc > 1) {
        len = stoull(argv[1]);
    }
    int_vector<> v(len, 7); // create a vector of length 50
    for (uint64_t i=1; i < v.size(); ++i) {
        uint64_t x = v[i-1];
        v[i] = (x*7-1) % 97;
    }
    cout << "size_in_bytes(v)=" << size_in_bytes(v) << endl;

    rmq_succinct_sct<> rmq(&v);
    cout << "size_in_bytes(rmq)/v.size()=" << ((double)size_in_bytes(rmq))/v.size() << endl;

    if (v.size() < 100) {   // if the vector is small
        cout << "v = " << v << endl; // output it
    }
    uint64_t left = 0, right = v.size()-1;

    // rmq does not need v to answer the queries
    util::clear(v); // so we can free the space for v

    uint64_t count = 0;

    while (left < right) {
        uint64_t min_pos = rmq(left, right);
        if (++count < 10) {
            cout << "minimum of v[" << left << "..." << right << "] at index ";
            cout << min_pos << endl;
        }
        if (min_pos - left > right - min_pos) {
            right = min_pos - 1;
        } else {
            left = min_pos + 1;
        }
    }
    cout << endl;
    cout << "write_structure<JSON_FORMAT>:" << endl;
    cout << "----------------------------- " << endl;
    write_structure<JSON_FORMAT>(rmq, cout);
    cout << endl;
}
