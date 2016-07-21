#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/bit_vectors.hpp>
#include <c++/4.9.3/bitset>

using namespace sdsl;
using namespace std;

int main()
{
    bit_vector b2 = bit_vector(10, 1);
    rank_support_v<1> b_rank(&b2); // <- pointer to b
    cout << "(" << 2 << ", " << b_rank(0) << ") ";

    bit_vector b = {0,1,0,1};
    rank_support_v<1> b_r1(&b);
    for (size_t i=0; i<=b.size(); ++i)
        cout << i << ": "<< b_r1(i) << endl;
}
