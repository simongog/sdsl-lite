#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    /*wt_huff_int<> wt;
    construct_im(wt, int_vector<>({1981, 1974, 1990, 1974, 2014, 1974}));
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << wt << endl;
    size_t idx = 5;
    auto r_c = wt.inverse_select(idx);
    cout << get<0>(r_c)+1 << " occurrence(s) of "
         << get<1>(r_c) << " in [0.." << idx << "]" << endl;

    for (size_t i=0; i < wt.size() and wt[i]!='\n'; ++i)
        cout << wt[i] << "\t";
    cout << endl;
    cout << "number of lines  : " << wt.rank(wt.size(), '\n') << endl;
    cout << "first '=' in line: " << wt.rank(wt.select(1, '='),'\n')+1 << endl;*/
}