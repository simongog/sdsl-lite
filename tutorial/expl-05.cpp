#include <iostream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    bit_vector b = bit_vector(80*(1<<20), 0);
    for (size_t i=0; i < b.size(); i+=100)
        b[i] = 1;
    cout << size_in_mega_bytes(b) << endl;
    rrr_vector<63> rrrb(b);
    cout << size_in_mega_bytes(rrrb) << endl;
    sd_vector<> sdb(b);
    cout << size_in_mega_bytes(sdb) << endl;
}
