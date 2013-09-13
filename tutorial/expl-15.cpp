#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <utility>

using namespace std;
using namespace sdsl;

template<class value_type, class size_type>
ostream& operator<<(ostream& os, const std::pair<value_type, size_type>& p)
{
    return os << "(" << p.first << "," << p.second << ")";
}

int main()
{
    wt_int<rrr_vector<63>> wt;
    construct_im(wt, "6   1000 1 4 7 3   18 6 3", 'd');
    auto res = wt.range_search_2d(1, 5, 4, 18);
    for (auto point : res.second)
        cout << point << " ";
    cout << endl;
}
