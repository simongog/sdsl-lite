#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    wt_huff_int<rrr_vector<63>> wt;
    construct_im(wt, int_vector<>({1981, 1974, 1990, 1974, 2014, 1974}));
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << wt << endl;
    size_t idx = 5;
    auto r_c = wt.inverse_select(idx);
    cout << get<0>(r_c)+1 << " occurrence(s) of "
         << get<1>(r_c) << " in [0.." << idx << "]" << endl;
}
