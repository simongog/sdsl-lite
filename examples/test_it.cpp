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

    std::vector<pair<uint64_t, uint64_t>> points = {{770,775},{39683,55811}};
    /*int size = 8;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++) {
            points.push_back(std::make_pair(i,j));
        }
    }*/

    std::sort(points.begin(), points.end(), sort_by_z_order());

    std::cout << "Sorted points" << std::endl;
    for (int i = 0; i < points.size(); ++i) {
        std::cout << points[i].first << "," << points[i].second << "\t";
    }

    std::cout << std::endl;
}
