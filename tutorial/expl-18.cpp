#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace sdsl;

int main()
{
    csa_wt<wt_huff<rrr_vector<63>>,4,8> csa; // 接尾辞配列
    construct(csa, "expl-18.cpp", 1);
    cout << "count(\"配列\") : " << count(csa, "配列") << endl;
    auto occs = locate(csa, "\n");
    sort(occs.begin(), occs.end());
    auto max_line_length = occs[0];
    for (size_t i=1; i < occs.size(); ++i)
        max_line_length = std::max(max_line_length, occs[i]-occs[i-1]+1);
    cout << "max line length : " << max_line_length << endl;
}
