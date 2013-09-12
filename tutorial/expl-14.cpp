#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    wt_blcd<rrr_vector<63>> wt;
    construct_im(wt, "ハローワールド！", 1);
    for (size_t i=0; i < wt.size(); ++i)
        cout << wt[i];
    cout << endl;
    wt_hutu<>::size_type cnt1 = 0, cnt2 = 0, _;
    wt.lex_count(0, wt.size(), 0x80, cnt1, _);
    wt.lex_count(0, wt.size(), 0xbf, _, cnt2);
    cout << "# of chars        : " << wt.size() << endl;
    cout << "# of UTF-8 symbols: " << cnt1 + cnt2 << endl;
}
