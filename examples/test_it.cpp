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

    //std::vector<pair<uint64_t, uint64_t>> points = {{770,775},{39683,55811}};

    std::vector<pair<uint, uint>> points;
    int size = 8;
    const int k = 2;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++) {
            points.push_back(std::make_pair(i,j));
        }
    }

    k2_treap<k, bit_vector> tree;
    tree.set_height(log(size)/log(k));

    std::sort(points.begin(), points.end(), [&](const std::pair<uint32_t, uint32_t>& lhs, const std::pair<uint32_t, uint32_t>& rhs){
        return tree.sort_by_z_order(lhs, rhs);
    });

    vector<pair<uint32_t, uint32_t>> correct_order = {{0,0},{0,1},{1,0},{1,1},{0,2},{0,3},{1,2},{1,3},
                                                      {2,0},{2,1},{3,0},{3,1},{2,2},{2,3},{3,2},{3,3},
                                                      {0,4},{0,5},{1,4},{1,5},{0,6},{0,7},{1,6},{1,7},
                                                      {2,4},{2,5},{3,4},{3,5},{2,6},{2,7},{3,6},{3,7},
                                                      {4,0},{4,1},{5,0},{5,1},{4,2},{4,3},{5,2},{5,3},
                                                      {6,0},{6,1},{7,0},{7,1},{6,2},{6,3},{7,2},{7,3},
                                                      {4,4},{4,5},{5,4},{5,5},{4,6},{4,7},{5,6},{5,7},
                                                      {6,4},{6,5},{7,4},{7,5},{6,6},{6,7},{7,6},{7,7},
    };

    /*int size = 8;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++) {
            points.push_back(std::make_pair(i,j));
        }
    }*/

    //std::sort(points.begin(), points.end(), sort_by_z_order());

    std::cout << "Sorted points" << std::endl;
    for (uint i = 0; i < points.size(); ++i) {
        std::cout << points[i].first << "," << points[i].second << "\t";
    }

    std::cout << std::endl;
}
