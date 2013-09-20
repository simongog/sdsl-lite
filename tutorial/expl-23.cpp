#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    cst_sct3<csa_wt<wt_huff<rrr_vector<>>>, lcp_support_sada<>> cst1;
    construct(cst1, "english.200MB", 1);
    cout << "cst1.lcp in MiB : " << size_in_mega_bytes(cst1.lcp) << endl;
    util::clear(cst1);
    cst_sct3<csa_wt<wt_huff<rrr_vector<>>>, lcp_dac<>> cst2;
    construct(cst2, "english.200MB", 1);
    cout << "cst2.lcp in MiB : " << size_in_mega_bytes(cst2.lcp) << endl;
}
