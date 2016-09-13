#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <bitset>

using namespace std;
using namespace sdsl;

int main()
{
/*
    wt_huff_int<> wt;
    int_vector<32> vec = {3,12,4,4,5,1,6,4,2};
    auto asd = vec.raw;

    std::cout << "Serialized vec incoming" << std::endl;
    store_to_file(vec.raw, "/home/d056848/myfile");
    std::cout << "Serialized vec incoming" << std::endl;
    construct_im(wt, vec.raw, 4);
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << wt << endl;


    for (size_t i=0; i < wt.size() and wt[i]!='\n'; ++i)
        cout << std::to_string(wt[i]) << "\t";
    cout << endl;
*/


    int_vector<64> x_vec(6, 0);
    util::set_to_id(x_vec);
    cout << x_vec << endl;  // 0 1 2 3 4 5
    wt_huff_int<> wt;
    construct_im(wt, x_vec.raw, 8);
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << "wt.size : " << wt.size() << endl;
    cout << wt << endl; // 0 1 2 3 4 5*/
}