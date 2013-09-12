#include <iostream>
#include <sdsl/vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    int_vector<> v = {3,2,1,0,2,1,3,4,1,1,1,3,2,3};
    v[1] = 0;
    util::bit_compress(v);
    cout << v << endl;
    cout << size_in_bytes(v) << endl;
}
