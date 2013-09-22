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
    cout << bps.find_close(0) << ", "
         << bps.find_open(3) << ", "
         << bps.enclose(4) << ", "
         << bps.double_enclose(13, 16) << endl;
}
