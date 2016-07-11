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

    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    void get_paper_k2_tree(k2_treap<t_k, t_bv, t_rank> &idx){
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
        std::string tmp_prefix = ram_file_name("k2_treap_test");
        std::vector<std::pair<uint32_t, uint32_t>> coordinates = {{0,1},{1,2},{1,3},{1,4},{7,6},{8,6},{8,9},{9,6},{9,8},{9,10},{10,6},{10,9}};
        k2_treap<t_k, t_bv, t_rank> tmp(coordinates, tmp_prefix);
        tmp.swap(idx);
    }

    TEST(K2TreapInternalTest, testZOrderSort) {

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

    TEST(K2TreapInternalTest, testZOrder100) {

        std::vector<pair<uint, uint>> points = {{30286, 49131},
                                                {62359, 18315},
                                                {30605, 48934},
                                                {2973,  38167},
                                                {63409, 48393},
                                                {49407, 19780},
                                                {14211, 19800},
                                                {33583, 2142},
                                                {44435, 10393},
                                                {35375, 59222},
                                                {65095, 1657},
                                                {58958, 56987},
                                                {56664, 57987},
                                                {1076,  20759},
                                                {468,   47602},
                                                {46582, 27220},
                                                {43537, 42103},
                                                {63839, 60066},
                                                {48772, 60560},
                                                {34521, 56903},
                                                {52853, 30779},
                                                {16605, 42713},
                                                {23013, 45548},
                                                {5812,  30437},
                                                {45365, 40152},
                                                {45916, 31791},
                                                {26603, 64346},
                                                {39982, 6616},
                                                {8557,  1274},
                                                {11990, 38348},
                                                {2051,  42664},
                                                {63997, 25805},
                                                {5562,  54495},
                                                {5700,  53280},
                                                {37470, 60306},
                                                {49416, 51532},
                                                {8595,  25213},
                                                {35894, 30012},
                                                {34784, 53726},
                                                {35926, 36810},
                                                {9875,  9859},
                                                {18353, 64753},
                                                {36388, 53738},
                                                {12072, 51156},
                                                {15945, 35438},
                                                {27097, 42045},
                                                {37067, 40719},
                                                {41221, 64392},
                                                {19933, 15127},
                                                {14572, 46869},
                                                {15133, 33072},
                                                {64546, 60344},
                                                {58441, 30066},
                                                {19512, 54802},
                                                {18908, 1706},
                                                {10694, 16841},
                                                {5820,  28620},
                                                {4337,  52281},
                                                {12870, 31461},
                                                {63761, 40465},
                                                {49393, 1154},
                                                {22513, 12605},
                                                {10068, 39068},
                                                {27097, 14940},
                                                {48617, 59872},
                                                {9458,  36990},
                                                {54599, 16561},
                                                {9385,  16459},
                                                {40688, 4998},
                                                {1748,  46513},
                                                {52286, 59351},
                                                {62820, 56632},
                                                {7832,  42385},
                                                {44134, 47311},
                                                {24769, 13050},
                                                {12919, 23618},
                                                {26120, 31921},
                                                {41343, 3623},
                                                {28967, 23315},
                                                {14419, 30449},
                                                {38136, 18051},
                                                {17793, 58364},
                                                {29281, 36824},
                                                {23175, 33272},
                                                {64569, 54760},
                                                {17632, 5750},
                                                {27310, 20531},
                                                {32301, 52174},
                                                {38742, 44580},
                                                {59689, 10151},
                                                {53465, 1210},
                                                {33412, 3113},
                                                {42167, 44717},
                                                {60541, 56972},
                                                {64630, 35365},
                                                {17528, 38430},
                                                {45583, 51243},
                                                {35222, 7964},
                                                {34983, 3252},
                                                {23070, 7096}};
        int size = 8;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                points.push_back(std::make_pair(i, j));
            }
        }

        std::sort(points.begin(), points.end(), sort_by_z_order());


        std::cout << "#########Links##############"<< std::endl;
        for (int n = 0; n < points.size(); ++n) {
            std::cout << "{"<< points[n].first << "," << points[n].second << "},";
        }
        std::cout << std::endl;
    }

    TEST(K2TreapInternalTest, test_calculate_subtree_number_and_new_relative_coordinates) {
        k2_treap<2, bit_vector> tree;
        get_paper_k2_tree(tree);
        //crappy test with magic numbers ;-)
        std::pair<uint, uint> test_link = std::make_pair(8,9);
        uint subtree_number = tree.calculate_subtree_number_and_new_relative_coordinates(test_link,0);
        ASSERT_EQ(subtree_number,3);
        ASSERT_EQ(test_link.first,0);
        ASSERT_EQ(test_link.second,1);
        uint subtree_number_2 = tree.calculate_subtree_number_and_new_relative_coordinates(test_link,1);
        ASSERT_EQ(subtree_number_2,0);
    }
}
// namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
