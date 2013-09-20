#include <iostream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    bit_vector b = {1,1,0,1,0,0,1};
    cout << b << endl;
    b = bit_vector(80*(1<<20), 0);
    for (size_t i=0; i < b.size(); i+=100)
        b[i] = 1;
    cout << size_in_mega_bytes(b) << endl;
}
