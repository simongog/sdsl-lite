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
        ofstream out("data.sdsl");
        uint64_t x = 42;
        serialize(x, out);
    }
    {
        ifstream in("data.sdsl");
        uint64_t x = 0;
        load(x, in);
        cout << x << endl;
    }
}
