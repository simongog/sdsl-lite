#include <iostream>
#include <sdsl/vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    int_vector<> v(10000000, 3);
    v[0] = 1ULL<<63;
    util::bit_compress(v);
    cout << size_in_mega_bytes(v) << endl;
    vlc_vector<> vv(v);
    cout << size_in_mega_bytes(vv) << endl;
}
