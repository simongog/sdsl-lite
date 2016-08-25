#include <iostream>
#include <vector>
#include <tuple>
#include <complex>
#include <sdsl/k2_tree_algorithm.hpp>

using namespace sdsl;
using namespace std;

void print(std::vector<uint32_t>& to_print, uint32_t source_node){
    std::cout << "Reachable nodes from " << source_node <<std::endl;
    for (auto& link: to_print){
        std::cout << link << "\t";
    }
    std::cout << std::endl;
}

int main()
{
    typedef k2_tree<2, bit_vector, bit_vector, false, 2> k2;
    typedef k2_tree_partitioned<2, k2, true> k2_rrr;


    // Initialize treap with a vector of (x,y,weight) elements
    vector<pair<uint32_t, uint32_t>> coords = {{0,0},{0,1},{1,2},{1,3},{1,4},{7,6},{8,6},{8,9},{9,6},{9,8},{9,10},{10,6},{10,9}};//{{0,0},{0,1},{1,4},{1,3},{7,6},{1,2}};
    k2_rrr k2treap;
    construct_im(k2treap, coords,10);

    cout << "Points in the k2treap: " << k2treap.size() << endl;

    std::vector<uint32_t> result;
/*
    if (k2treap.check_link(std::make_pair((uint)0,(uint)0))){
        std::cout << "1" << std::endl;
    } else {
        std::cout << "0" << std::endl;
    }
    */
    //k2treap.check_link_shortcut(std::make_pair((uint) 0,(uint) 5));

/*
    for (uint i = 0; i < 15; ++i) {
        for (uint j = 0; j < 15; ++j) {
            if (k2treap.check_link(std::make_pair(i,j))){
                std::cout << 1 << " ";
            } else {
                std::cout << 0 << " ";
            }
        }
        std::cout << "\r\n";
    }
*/
    //k2treap.direct_links((uint32_t) 4, result);
   // k2treap.direct_links_shortcut((uint)7, result);
    k2treap.direct_links2((uint)6, result);


    for (uint32_t i = 0; i < 16; i++){
        k2treap.direct_links2(i, result);
        print(result, i);
    }
/*
    std::cout << "inverse links" << std::endl;

    for (uint32_t i = 0; i < 16; i++){
        k2treap.inverse_links_shortcut(i, result);
        print(result, i);
    }*/

    /*store_to_file(k2treap, "tmp.tmp");

    k2_tree_hybrid<4,2,2,2, bit_vector, rrr_vector<63>> k2treap2;
    load_from_file(k2treap2, "tmp.tmp");

    sdsl::remove("tmp.tmp");*/

    /*auto range_it = range_3d(k2treap, {2,2}, {10,10});
    while (range_it) {
        auto p = *range_it;
        cout << real(p) <<" weight: "<< imag(p) << endl;
        ++range_it;
    }

    cout<<"---"<<endl;
    {
        k2_rrr k2t;
        vector<pair<uint32_t, uint32_t>> coords = {{1,2},{2,3},{3,1}};
        construct_im(k2t, coords);
        auto topk_it = range_3d(k2t, {0,0}, {10,10});
        while (topk_it) {
            auto p = *topk_it;
            cout << real(p) <<" weight: "<< imag(p) << endl;
            ++topk_it;
        }
    }*/
}
