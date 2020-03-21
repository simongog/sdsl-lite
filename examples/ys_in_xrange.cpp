#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>
#include <iostream>
#include <stack>

using namespace sdsl;
using namespace std;


int main()
{
    typedef wt_int<> t_wt;
    t_wt wt;
    construct_im(wt, "9 4 3 2 1 4 6 3 1 4 6 5 3 2 1 3 5 3 2 3 4",'d');
    cout << wt << endl;
    auto y_it = ys_in_x_range(wt, 0, wt.size());
    while (y_it) {
        auto y = *y_it;
        cout << get<0>(y) << " ("<< get<1>(y) << "," << get<2>(y) << ")" << endl;
        ++y_it;
    }

    cout << "count[0,"<<wt.size()-1<<"][1,16] = " << count(wt, {0, wt.size()-1}, {1,16}) << endl;
    cout << "count[0,"<<wt.size()-1<<"][1,8] = " << count(wt, {0, wt.size()-1}, {1,8}) << endl;
    cout << "count[0,"<<wt.size()-1<<"][2,3] = " << count(wt, {0, wt.size()-1}, {2,3}) << endl;


    cout << "map_to_sorted( {"<<wt<<"}, [5,13), [1,5) )" << endl;
    auto mts_it = map_to_sorted_sequence(wt, {5, 13}, {0,3});
    while (mts_it) {
        cout << get<0>(*mts_it);
        cout << " ["<<get<0>(get<1>(*mts_it)) <<","<<get<1>(get<1>(*mts_it))<<"]"<<endl;
        ++mts_it;
    }
}

