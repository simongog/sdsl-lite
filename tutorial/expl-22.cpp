#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    cst_sct3<csa_wt<wt_int<rrr_vector<>>>> cst;
    int_vector<> data(100000, 0, 10);
    for (size_t i=0; i < data.size(); ++i)
        data[i] = 1 + rand()%1023;
    construct_im(cst, data);
    cout << "cst.csa.sigma: " << cst.csa.sigma << endl;
    for (size_t k=0; k<3; ++k)
        cout << "H" << k << "(data) : " <<  get<0>(Hk(cst, k)) << endl;
}
