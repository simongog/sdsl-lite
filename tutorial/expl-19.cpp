#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace sdsl;

int main()
{
    csa_wt<wt_int<rrr_vector<63>>> csa;
    construct_im(csa, "1 8 15 23 1 8 23 11 8", 'd');
    cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", csa);
}
