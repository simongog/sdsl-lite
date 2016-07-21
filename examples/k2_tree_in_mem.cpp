#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree.hpp>
#include <sdsl/bit_vectors.hpp>

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
    typedef k2_tree<2,rrr_vector<63>> k2_rrr;
    k2_rrr k2treap;

    // Initialize treap with a vector of (x,y,weight) elements
    vector<pair<uint32_t, uint32_t>> coords = {{0,0},{0,1},{1,2},{1,3},{1,4},{7,6},{8,6},{8,9},{9,6},{9,8},{9,10},{10,6},{10,9}};//{{0,0},{0,1},{1,4},{1,3},{7,6},{1,2}};
    construct_im_bottom_up(k2treap, coords);

    cout << "Points in the k2treap: " << k2treap.size() << endl;

    cout << "Points in the rectangle from (2,1) to (3,3): " ;
    //cout << count(k2treap, {2,1}, {3,3}) << endl;

    cout << "Report all points in rectangle from (2,2) to (10,10)" << endl;
    std::vector<uint32_t> result;
/*
    k2treap.check_link(std::make_pair(9,6));


    std::pair<uint64_t, uint64_t> asd = std::make_pair(4,5);
    bool asd2 = k2treap.check_link(asd);

    for (int i = 0; i < 15; ++i) {
        for (int j = 0; j < 15; ++j) {
            if (k2treap.check_link(std::make_pair(i,j))){
                std::cout << 1 << " ";
            } else {
                std::cout << 0 << " ";
            }
        }
        std::cout << "\r\n";
    }

*/
    for (uint32_t i = 0; i < 16; i++){
        k2treap.direct_links2(i, result);
        print(result, i);
    }

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
