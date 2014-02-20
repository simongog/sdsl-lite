#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    wt_int<> wt;
    //0 1 2 3 4 5 6 7 8 9 0 1 2 3 4  5 6 7 8 9 0  1 2
    construct_im(wt, int_vector<> {1,2,2,4,5,3,2,4,5,3,2,4,7,4,2,52,7,4,2,1,5,74,3});

    cout << "wt = " << wt << endl;
    cout << "wt[0,3] intersect with wt[7,10]" << endl;
    auto res = intersect(wt, {{0,3},{7,10}});
    cout << "element : sum of occurrences in ranges" << endl;
for (auto x : res) {
        cout << x.first << " : " << x.second <<endl;
    }

}
