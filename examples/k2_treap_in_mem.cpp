#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_treap.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace sdsl;
using namespace std;

int main()
{
    typedef k2_treap<2,rrr_vector<63>> k2_rrr;
    k2_rrr k2treap;

    // Initialize treap with a vector of (x,y,weight) elements
    vector<pair<uint32_t, uint32_t>> coords = {{0,1},{1,4},{1,3},{7,6},{1,2}};
    construct_im(k2treap, coords);

    cout << "Points in the k2treap: " << k2treap.size() << endl;

    cout << "Points in the rectangle from (2,1) to (3,3): " ;
    cout << count(k2treap, {2,1}, {3,3}) << endl;

    cout << "Report all points in rectangle from (2,2) to (10,10)" << endl;
    auto range_it = range_3d(k2treap, {2,2}, {10,10});
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
    }
}
