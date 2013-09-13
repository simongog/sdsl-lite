#include <iostream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    bit_vector b = bit_vector(8000, 0);
    for (size_t i=0; i < b.size(); i+=100)
        b[i] = 1;
    rank_support_v<1> b_rank(&b);
    for (size_t i=0; i<=b.size(); i+= b.size()/4)
        cout << "(" << i << ", " << b_rank(i) << ") ";
    cout << endl;
}
