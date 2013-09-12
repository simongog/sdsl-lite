#include <iostream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    bit_vector b = {0,1,0,1,1,1,0,0,0,1,1};
    size_t cnt10 = rank_support_v<10,2>(&b)(b.size());
    select_support_mcl<10,2> b_sel10(&b);

    for (size_t i=1; i <= cnt10; ++i)
        cout << b_sel10(i) << " ";
    cout << endl;
}
