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
    typedef k2_treap<3,rrr_vector<63>> k2_rrr;
    k2_rrr k2treap;

    // Initialize treap with a vector of (x,y,weight) elements
    construct_im(k2treap, {{1,2,3},{2,2,6},{4,4,1},{3,3,2},{3,1,8}});

    cout << "Points in the k2treap: " << k2treap.size() << endl;

    cout << "Points in the rectangle from (2,1) to (3,3): " ;
    cout << count(k2treap, {2,1}, {3,3}) << endl;

    cout << "Heaviest points in rectangle from (0,0) to (2,8):" << endl;
    auto topk_it = top_k(k2treap, {0,0}, {2,8});
    while (topk_it) {
        auto point_weight = *topk_it;
        cout << point_weight.first <<" weight: "<<point_weight.second << endl;
        ++topk_it;
    }

    cout << "Report all points in rectangle from (2,2) to (10,10)" << endl;
    cout << "with weight in [2..6]:" << endl;
    auto range_it = range_3d(k2treap, {2,2}, {10,10}, {2,100});
    while (range_it) {
        auto point_weight = *range_it;
        cout << point_weight.first <<" weight: "<<point_weight.second << endl;
        ++range_it;
    }

    cout<<"---"<<endl;
    {
        k2_rrr k2t;
        construct_im(k2t, {{1,2,3},{2,3,3},{3,1,3}});
        auto topk_it = top_k(k2t, {0,0}, {10,10});
        while (topk_it) {
            auto point_weight = *topk_it;
            cout << point_weight.first <<" weight: "<<point_weight.second << endl;
            ++topk_it;
        }
    }
}
