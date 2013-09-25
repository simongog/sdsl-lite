#include <sdsl/wavelet_trees.hpp>
#include <sdsl/bit_vectors.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

template<class t_wt>
void show_size(string input, string description)
{
    t_wt wt;
    construct(wt, input, 1);
    cout << size_in_bytes(wt) << " - " << description << endl;
}

int main(int argc, char* argv[])
{
    show_size<wt_blcd<>>(argv[0], "wt_bcld with default arguments");
    show_size<wt_blcd<bit_vector_il<>>>(argv[0], "use bit_vector_il");
    show_size<wt_blcd<rrr_vector<>>>(argv[0], "use rrr_vector");
    show_size<wt_blcd<rrr_vector<63>>>(argv[0], "use rrr_vector<63>");
    show_size<wt_blcd<bit_vector>>(argv[0], "use bit_vector - default bitvector");
    show_size<wt_blcd<bit_vector, rank_support_v<>>>(argv[0], "use rank_support_v - default rank");
    show_size<wt_blcd<bit_vector, rank_support_v5<>>>(argv[0], "use rank_support_v5 - more space efficient");
    show_size<wt_blcd<bit_vector, rank_support_scan<>>>(argv[0], "use rank_support_scan - no overhead for rank");
    show_size<wt_blcd<bit_vector, rank_support_v<>, select_support_mcl<1>>>(argv[0], "use select_support_mcl<1> - default select 1");
    show_size<wt_blcd<bit_vector, rank_support_v<>, select_support_scan<1>>>(argv[0], "use select_support_scan - no overhead for select 1");
    show_size<wt_blcd<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>>>(argv[0], "use select_support_mcl<0> - default select 0");
    show_size<wt_blcd<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_scan<0>>>(argv[0], "use select_support_scan - no overhead for select 0");
    return 0;
}

