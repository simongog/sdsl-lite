
#include <sdsl/bp_support.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    //              0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
    //              ( ( ( ) ( ) ) ( ) ( ( ( ) ( ) ) ( ) ) )
    bit_vector b = {1,1,1,0,1,0,0,1,0,1,1,1,0,1,0,0,1,0,0,0};
    bp_support_sada<> bps(&b); // <- pointer to b
    for (size_t i=0; i < b.size(); ++i)
        cout << bps.excess(i)<< " ";
    cout << endl;
    cout << bps.rank(0) << ", " // inclusive rank for BPS!!!
         << bps.select(4) << endl;
}
