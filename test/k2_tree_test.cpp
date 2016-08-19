#include "sdsl/k2_tree.hpp"
#include "gtest/gtest.h"
// #include <vector>
// #include <tuple>
// #include <string>
// #include <algorithm> // for std::min. std::sort
// #include <random>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;

template<class T>
class k2_tree_test_k_2 : public ::testing::Test { };

template<class T>
class k2_tree_test_k_5 : public ::testing::Test { };

using testing::Types;

namespace k2_tree_test_nm
{
    template<typename t_tree>
    void check_t_l(t_tree& tree, vector<unsigned> expected_t,
                   vector<unsigned> expected_l)
    {
    ASSERT_EQ(expected_t.size(), tree.get_t().size());
    ASSERT_EQ(expected_l.size(), tree.get_l().size());
    for(unsigned i = 0; i < expected_t.size(); i++)
        ASSERT_EQ(expected_t[i], tree.get_t().get_int(i, 1));
    for(unsigned i = 0; i < expected_l.size(); i++)
        ASSERT_EQ(expected_l[i], tree.get_l().get_int(i, 1));
    }
};

typedef Types<
    k2_tree<2, bit_vector, rank_support_v<>>> k_2_implementations;

TYPED_TEST_CASE(k2_tree_test_k_2, k_2_implementations);

TYPED_TEST(k2_tree_test_k_2, build_from_matrix_test)
{
    vector<vector <int>> mat({{1, 1, 0, 0},
                              {0, 1, 0, 0},
                              {0, 0, 1, 1},
                              {0, 0, 1, 0}});

    TypeParam tree(mat);
    vector<unsigned> expected_l = {1,1,0,1,1,1,1,0};
    k2_tree_test_nm::check_t_l(tree, {1, 0, 0 ,1}, expected_l);

    mat = vector<vector <int>> ({{0, 0, 0, 0},
                                 {0, 0, 0, 0},
                                 {0, 0, 0, 0},
                                 {0, 0, 0, 0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    mat = vector<vector<int>>({{0, 0},
                               {0, 0}});
    tree = TypeParam(mat);
    ASSERT_TRUE(tree.get_t().empty());
    ASSERT_TRUE(tree.get_l().empty());

    // Size is minor than k:
    mat = vector<vector<int>>({{0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});


    // Size is non a power of k:
    mat = vector<vector<int>>({{0, 0, 1},
                               {0, 1, 0},
                               {0, 1, 0}});
    tree = TypeParam(mat);
    expected_l = {0,0,0,1,1,0,0,0,0,1,0,0};
    k2_tree_test_nm::check_t_l(tree, {1, 1, 1 ,0}, expected_l);

    mat = vector<vector <int>>({{0, 0, 0},
                                {1, 0, 1},
                                {0, 1, 1}});
    tree = TypeParam(mat);
    expected_l = {0, 0, 1, 0,  0, 0, 1, 0,  0, 1, 0, 0,  1, 0, 0, 0};
    k2_tree_test_nm::check_t_l(tree, {1, 1, 1 ,1}, expected_l);

    // Sample from 'k^2 trees for compact web graph representation' paper
    mat = vector<vector<int>>({{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0},
                               {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1},
                               {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0}});
    tree = TypeParam(mat);
    vector<unsigned> expected_t = {1, 0, 1, 1,  1, 1, 0, 1,  0, 1, 0, 0,
                                   1, 0, 0, 0,  1, 1, 0, 0,  1, 0, 0, 0,
                                   0, 0, 0, 1,  0, 1, 0, 1,  1, 1, 1, 0};

    expected_l = {0, 1, 0, 0,  0, 0, 1, 1,  0, 0, 1, 0,  0, 0, 1, 0,
                  1, 0, 1, 0,  1, 0, 0, 0,  0, 1, 1, 0,  0, 0, 1, 0,
                  0, 1, 0, 0};
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);
}

TYPED_TEST(k2_tree_test_k_2, build_from_edges_array)
{
    vector<std::tuple<typename TypeParam::idx_type, typename TypeParam::idx_type>> e;
    // std::tuple<typename TypeParam::idx_type, typename TypeParam::idx_type> t({1, 2});
    e.push_back(std::make_tuple(typename TypeParam::idx_type(1),
                                typename TypeParam::idx_type(2))); // ({{1, 2}, {1, 3}});
    TypeParam tree(e, 4);
    // vector<unsigned> expected_l = {1,1,0,1,1,1,1,0};
    // k2_tree_test_nm::check_t_l(tree, {1, 0, 0 ,1}, expected_l);

    // mat = vector<vector <int>> ({{0, 0, 0, 0},
                                 // {0, 0, 0, 0},
}

TYPED_TEST(k2_tree_test_k_2, neighbors_test)
{
    vector<vector <int>> mat({{1, 1, 0, 0},
                              {0, 1, 0, 0},
                              {0, 0, 1, 1},
                              {0, 0, 1, 0}});

    TypeParam tree(mat);
    auto neigh_0 = tree.neigh(0);
    vector<unsigned>expected_neigh_0({0, 1});
    ASSERT_EQ(expected_neigh_0.size(), neigh_0.size());
    for(unsigned i = 0; i < neigh_0.size(); i++)
        ASSERT_EQ(expected_neigh_0[i], neigh_0[i]);

    auto neigh_3 = tree.neigh(3);
    vector<unsigned>expected_neigh_3({2});
    ASSERT_EQ(expected_neigh_3.size(), neigh_3.size());
    for(unsigned i = 0; i < neigh_3.size(); i++)
        ASSERT_EQ(expected_neigh_3[i], neigh_3[i]);

    mat = vector<vector <int>>({{0, 0, 0},
                                {1, 0, 1},
                                {0, 1, 1}});
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0.size());

    auto neigh_1 = tree.neigh(1);
    auto expected_neigh_1 = vector<unsigned>({0, 2});
    ASSERT_EQ(expected_neigh_1.size(), neigh_1.size());
    for(unsigned i = 0; i < neigh_1.size(); i++)
        ASSERT_EQ(expected_neigh_1[i], neigh_1[i]);


    mat = vector<vector <int>>({{0, 0},
                                {0, 0}});
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0.size());
}

TYPED_TEST(k2_tree_test_k_2, reverse_neighbors_test)
{
    vector<vector <int>> mat({{1, 0, 0, 0, 1},
                              {0, 0, 0, 0, 0},
                              {0, 0, 1, 1, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 1, 0, 1}});

    auto tree = TypeParam(mat);
    auto r_neigh_0 = tree.reverse_neigh(0);
    auto expected_r_neigh_0 = vector<unsigned>({0});
    auto r_neigh_1 = tree.reverse_neigh(1);
    auto r_neigh_2 = tree.reverse_neigh(2);
    auto expected_r_neigh_2 = vector<unsigned>({2, 4});
    ASSERT_EQ(expected_r_neigh_0.size(), r_neigh_0.size());
    ASSERT_EQ(0u, r_neigh_1.size());
    ASSERT_EQ(expected_r_neigh_2.size(), r_neigh_2.size());

    for(unsigned i = 0; i < r_neigh_0.size(); i++)
        ASSERT_EQ(expected_r_neigh_0[i], r_neigh_0[i]);

    for(unsigned i = 0; i < r_neigh_2.size(); i++)
        ASSERT_EQ(expected_r_neigh_2[i], r_neigh_2[i]);

    mat = vector<vector <int>>({{0, 0},
                                {0, 0}});
    tree = TypeParam(mat);
    r_neigh_0 = tree.reverse_neigh(0);
    ASSERT_EQ(0u, r_neigh_0.size());
}

TYPED_TEST(k2_tree_test_k_2, adj_test)
{
    vector<vector <int>> mat({{1, 0, 0, 0, 1},
                              {0, 0, 0, 0, 0},
                              {0, 0, 1, 1, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 1, 0, 1}});

    auto tree = TypeParam(mat);
    ASSERT_TRUE(tree.adj(0, 0));
    ASSERT_TRUE(tree.adj(0, 4));
    ASSERT_FALSE(tree.adj(4, 0));
    ASSERT_FALSE(tree.adj(7, 7));
    ASSERT_FALSE(tree.adj(1, 1));
    ASSERT_TRUE(tree.adj(2, 2));
    ASSERT_TRUE(tree.adj(2, 3));

    mat = vector<vector <int>>({{0}});
    tree = TypeParam(mat);
    ASSERT_FALSE(tree.adj(0,0));
    mat = vector<vector <int>>({{1}});
    tree = TypeParam(mat);
    ASSERT_TRUE(tree.adj(0,0));
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
