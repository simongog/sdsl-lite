#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <bitset>

using namespace std;
using namespace sdsl;

int main()
{
    wt_huff_int<> wt;
    int_vector<> vec = int_vector<>(9, 0, 64);
    util::set_to_id(vec);  // 0 1 2 3 4 5 6 7 8
    //{3,12,4,4,5,1,6,4,2};
    construct_im(wt, vec, 8);
    cout << wt << endl;

/*
    wt_huff_int<> wt;
    int_vector<> vec = int_vector<>(9, 0);
    util::set_to_id(vec);  // 0 1 2 3 4 5 6 7 8
    util::bit_compress(vec);
    //{3,12,4,4,5,1,6,4,2};
    construct_im(wt, vec);
    cout << wt << endl;


    for (size_t i=0; i < wt.size() and wt[i]!='\n'; ++i)
        cout << std::to_string(wt[i]) << "\t";
    cout << endl;
*/

/*
    int_vector<64> x_vec(6, 0);
    util::set_to_id(x_vec);
    cout << x_vec << endl;  // 0 1 2 3 4 5
    wt_huff_int<> wt;
    construct_im(wt, x_vec.raw, 8);
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << "wt.size : " << wt.size() << endl;
    cout << wt << endl; // 0 1 2 3 4 5*/

}