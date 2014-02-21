#include <sdsl/vectors.hpp>
#include <fstream>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    cout << (has_serialize<bit_vector>::value)  << endl;
    cout << (has_serialize<uint64_t>::value)  << endl;

    {
        ofstream out("data.sdsl");
        bit_vector b = {1,0,0,0,1,1,0};
        serialize(b, out);
    }
    {
        ifstream in("data.sdsl");
        bit_vector b;
        load(b, in);
        cout << b << endl;
    }
    {
        uint64_t x = 42;
        store_to_file(x, "data.sdsl");
    }
    {
        uint64_t x = 0;
        load_from_file(x, "data.sdsl");
        cout << x << endl;
    }
    {
        std::vector<uint32_t> x(10, 5);
        store_to_file(x, "data.sdsl");
    }
    {
        std::vector<uint32_t> x;
        load_from_file(x, "data.sdsl");
        cout << x.size() << endl;
        for (size_t i=0; i<x.size(); ++i)
            cout << x[i] << endl;
    }
}
