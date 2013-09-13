#include <iostream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    bit_vector b = {0,1,0,1};
    rank_support_v<1> b_r1(&b);
    rank_support_v<0> b_r0(&b);
    rank_support_v<10,2> b_r10(&b);
    rank_support_v<01,2> b_r01(&b);
    for (size_t i=0; i<=b.size(); ++i)
        cout << i << ": "<< b_r1(i) << " " << b_r0(i)
             << " " << b_r10(i) << " " << b_r01(i) << endl;
}
