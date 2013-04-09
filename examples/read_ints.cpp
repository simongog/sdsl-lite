#include <sdsl/suffix_trees.hpp>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;
using namespace sdsl;

int main()
{
    {
        typedef csa_wt<wt_int<>, 32, 64, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> > t_csa;
        typedef cst_sct3<t_csa> t_cst;
        t_cst cst;
        std::streamsize BUF_SIZE=4096;
        char* text = new char[BUF_SIZE];
        while (cin.getline(text, BUF_SIZE)) {
            construct_im(cst, (const char*)text, 'd');
            cout << "   i  SA LCP     T" << endl;
            for (uint64_t i=0; i <cst.csa.size(); ++i) {
                cout << std::setw(4) << i << std::setw(3) << cst.csa[i] << std::setw(3) << cst.lcp[i];
                for (uint64_t j=cst.csa[i]; j < cst.csa.size(); ++j) {
                    cout << (uint64_t)cst.csa.text[j];
                }
                cout << endl;
            }
        }
    }
}
