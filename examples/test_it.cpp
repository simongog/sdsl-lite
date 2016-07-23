#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/bit_vectors.hpp>
#include <c++/4.9.3/bitset>
#include <sdsl/k2_tree.hpp>

using namespace sdsl;
using namespace std;

template<uint8_t t_k,
        typename t_bv,
        typename t_rank>
void get_paper_k2_tree(k2_tree<t_k, t_bv, t_rank> &idx, uint access_shortcut_size = 0){
    //This Matrix is the same as the 16x16 matrix in the paper, the rows containing only 0s are omitted
    /*
    * 0 1 0 0 | 0 0 0 0 | 0 0 0
    * 0 0 1 1 | 1 0 0 0 | 0 0 0
    * 0 0 0 0 | 0 0 0 0 | 0 0 0
    * 0 0 0 0 | 0 0 0 0 | 0 0 0
    * -------------------------
    * 0 0 0 0 | 0 0 0 0 | 0 0 0
    * 0 0 0 0 | 0 0 0 0 | 0 0 0
    * 0 0 0 0 | 0 0 0 0 | 0 0 0
    * 0 0 0 0 | 0 0 1 0 | 0 0 0
    * -------------------------
    * 0 0 0 0 | 0 0 1 0 | 0 1 0
    * 0 0 0 0 | 0 0 1 0 | 1 0 1
    * 0 0 0 0 | 0 0 1 0 | 0 1 0
    */
    std::string tmp_prefix = ram_file_name("k2_tree_test");
    std::vector<std::pair<uint32_t, uint32_t>> coordinates = {{0,1},{1,2},{1,3},{1,4},{7,6},{8,6},{8,9},{9,6},{9,8},{9,10},{10,6},{10,9}};
    k2_tree<t_k, t_bv, t_rank> tmp(coordinates, tmp_prefix, false, access_shortcut_size);
    tmp.swap(idx);
}

int main()
{
    using namespace k2_treap_ns;
    uint access_shortcut_size = 3;
    k2_tree<2, bit_vector> tree;
    get_paper_k2_tree(tree, access_shortcut_size);

    //crappy test with magic numbers ;-)
    node_type* node = tree.check_link_shortcut((uint32_t)9, (uint32_t) 6);

    if (node != nullptr){
        std::cout << std::to_string(node->p.real())  << "sould be 8 \n";
        std::cout << std::to_string(node->p.imag())  << "sould be 6 \n";
        std::cout << std::to_string(node->idx)  << "sould be 52 \n";
    } else {
        std::cout << "No such link" << std::endl;
    }

}
