#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    wt_hutu<rrr_vector<63>> wt;
    construct_im(wt, "こんにちは世界", 1);
    for (size_t i=0; i < wt.size(); ++i)
        cout << wt[i];
    cout << endl;
    auto t1 = wt.lex_count(0, wt.size(), 0x80);
    auto t2 = wt.lex_count(0, wt.size(), 0xbf);
    cout << "# of bytes        : " << wt.size() << endl;
    cout << "# of UTF-8 symbols: " << get<1>(t1) + get<2>(t2) << endl;
}
