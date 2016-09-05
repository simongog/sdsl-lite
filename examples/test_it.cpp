#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <bitset>

using namespace std;
using namespace sdsl;

int main()
{
    bit_vector vec = {0,1,1};
    std::cout << std::bitset<8>(vec.get_int(0,8)) << std::endl;
    /*std::cout << std::bitset<8>(vec.get_int(0,1)) << std::endl;
    std::cout << std::bitset<8>(vec.get_int(1,2)) << std::endl;
    std::cout << std::bitset<8>(vec.get_int(1,5)) << std::endl;

    /*
    wt_huff_int<> wt;
    int_vector<8> vec = {3,12,4,4,5,1,6,4,2};
    construct_im(wt, vec);
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << wt << endl;

    /*
    for (size_t i=0; i < wt.size() and wt[i]!='\n'; ++i)
        cout << std::to_string(wt[i]) << "\t";
    cout << endl;

    int_vector<64> x_vec(6, 0);
    util::set_to_id(x_vec);
    cout << x_vec << endl;  // 0 1 2 3 4 5
    wt_huff_int<rrr_vector<63>> wt;
    construct_im(wt, x_vec.raw, 8);
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << "wt.size : " << wt.size() << endl;
    cout << wt << endl; // 0 1 2 3 4 5*/
}