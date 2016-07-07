#include "sdsl/k2_treap.hpp"
#include "sdsl/bit_vectors.hpp"
#include "gtest/gtest.h"
#include "../../../../../usr/local/include/c++/4.9.3/iterator"
#include <vector>
#include <tuple>
#include <string>
#include <algorithm> // for std::min. std::sort
#include <random>

namespace {

    using namespace sdsl;
    using namespace std;

    typedef int_vector<>::size_type size_type;

    template<class T>
    class k2_treap_test : public ::testing::Test {
    };

    TEST(K2TreapInternalTest, testZOrderSort) {
        k2_treap<3, bit_vector> treap;

        std::vector<pair<uint, uint>> points;
        int size = 8;
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++) {
                points.push_back(std::make_pair(i,j));
            }
        }

        std::sort(points.begin(), points.end(), sort_by_z_order());

        vector<pair<uint32_t, uint32_t>> correct_order = {{0,0},{0,1},{1,0},{1,1},{0,2},{0,3},{1,2},{1,3},
                                                          {2,0},{2,1},{3,0},{3,1},{2,2},{2,3},{3,2},{3,3},
                                                          {0,4},{0,5},{1,4},{1,5},{0,6},{0,7},{1,6},{1,7},
                                                          {2,4},{2,5},{3,4},{3,5},{2,6},{2,7},{3,6},{3,7},
                                                          {4,0},{4,1},{5,0},{5,1},{4,2},{4,3},{5,2},{5,3},
                                                          {6,0},{6,1},{7,0},{7,1},{6,2},{6,3},{7,2},{7,3},
                                                          {4,4},{4,5},{5,4},{5,5},{4,6},{4,7},{5,6},{5,7},
                                                          {6,4},{6,5},{7,4},{7,5},{6,6},{6,7},{7,6},{7,7},
        };

        ASSERT_EQ(points.size(),correct_order.size());
        std::cout << std::endl;
        for (int k = 0; k < correct_order.size(); ++k) {
            ASSERT_EQ(correct_order[k], points[k]);
        }
    }
}
// namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
