#include "sdsl/k2_tree.hpp"
#include "gtest/gtest.h"

#include <sstream>
#include <tuple>
#include <vector>

namespace
{

using namespace sdsl;
using namespace std;
using namespace k2_tree_ns;

typedef int_vector<>::size_type size_type;

template <class T>
class k2_tree_test_k_2 : public ::testing::Test{};

template <class T>
class k2_tree_test_k_3 : public ::testing::Test{};

template <class T>
class k2_tree_test : public ::testing::Test{};

template <class T>
class k2_tree_test_marked : public ::testing::Test{};

using testing::Types;

namespace k2_tree_test_nm
{
template <typename t_tree>
void check_t_l(t_tree &tree, vector<unsigned> expected_t,
               vector<unsigned> expected_l)
{
    ASSERT_EQ(expected_t.size(), tree.t().size());
    ASSERT_EQ(expected_l.size(), tree.l().size());
    for (unsigned i = 0; i < expected_t.size(); i++)
        ASSERT_EQ(expected_t[i], tree.t().get_int(i, 1));
    for (unsigned i = 0; i < expected_l.size(); i++)
        ASSERT_EQ(expected_l[i], tree.l().get_int(i, 1));
}

template <typename t_tree>
void assert_eq_tree(t_tree &tree1, t_tree &tree2)
{
    ASSERT_EQ(tree1.t().size(), tree2.t().size());
    ASSERT_EQ(tree1.l().size(), tree2.l().size());
    for (unsigned i = 0; i < tree1.t().size(); i++)
        ASSERT_EQ(tree1.t().get_int(i, 1), tree2.t().get_int(i, 1));

    for (unsigned i = 0; i < tree1.l().size(); i++)
        ASSERT_EQ(tree1.l().get_int(i, 1), tree2.l().get_int(i, 1));
}

template <typename t_tree>
void check_serialize_load(t_tree &tree)
{
    auto unserialized_tree = t_tree();
    std::stringstream ss;
    tree.serialize(ss);
    unserialized_tree.load(ss);
    ASSERT_EQ(tree, unserialized_tree);
}
}; // namespace k2_tree_test_nm

typedef Types<
    k2_tree<2, bit_vector, rank_support_v<>>,
    k2_tree<2, bit_vector>>
    k_2_implementations;

typedef Types<
    k2_tree<3, bit_vector, rank_support_v<>>,
    k2_tree<3, bit_vector>>
    k_3_implementations;

typedef Types<
    k2_tree<2, bit_vector>,
    k2_tree<3, bit_vector>,
    k2_tree<7, bit_vector>,
    k2_tree<2, rrr_vector<63>>,
    k2_tree<3, rrr_vector<63>>,
    k2_tree<5, bit_vector, rank_support_v<>>,
    k2_tree<4, bit_vector, rank_support_v<>>>
    Implementations;

typedef Types<
    k2_tree<2, bit_vector>,
    k2_tree<3, bit_vector>,
    k2_tree<7, bit_vector>,
    k2_tree<5, bit_vector, rank_support_v<>>,
    k2_tree<4, bit_vector, rank_support_v<>>>
    Implementations_bit_vector;

TYPED_TEST_CASE(k2_tree_test_k_2, k_2_implementations);

TYPED_TEST(k2_tree_test_k_2, build_from_matrix_test)
{
    vector<vector<int>> mat({{1, 1, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 1, 1},
                             {0, 0, 1, 0}});

    TypeParam tree(mat);
    vector<unsigned> expected_l = {1, 1, 0, 1, 1, 1, 1, 0};
    k2_tree_test_nm::check_t_l(tree, {1, 0, 0, 1}, expected_l);

    mat = vector<vector<int>>({{0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    mat = vector<vector<int>>({{0, 0},
                               {0, 0}});
    tree = TypeParam(mat);
    ASSERT_TRUE(tree.t().empty());
    ASSERT_TRUE(tree.l().empty());

    // Size is minor than k:
    mat = vector<vector<int>>({{0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    mat = vector<vector<int>>({{1}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0});

    // Size is non a power of k:
    mat = vector<vector<int>>({{0, 0, 1},
                               {0, 1, 0},
                               {0, 1, 0}});
    tree = TypeParam(mat);
    expected_l = {0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0};

    ASSERT_EQ(tree.get_number_edges(), 3);
    k2_tree_test_nm::check_t_l(tree, {1, 1, 1, 0}, expected_l);

    mat = vector<vector<int>>({{0, 0, 0},
                               {1, 0, 1},
                               {0, 1, 1}});
    tree = TypeParam(mat);
    expected_l = {0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0};
    k2_tree_test_nm::check_t_l(tree, {1, 1, 1, 1}, expected_l);

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
    vector<unsigned> expected_t = {1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0,
                                   1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0,
                                   0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0};

    expected_l = {0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0,
                  1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0,
                  0, 1, 0, 0};
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);
}

TYPED_TEST(k2_tree_test_k_2, build_from_edges_array)
{
    typedef std::tuple<typename TypeParam::idx_type,
                       typename TypeParam::idx_type>
        t_tuple;
    vector<std::tuple<typename TypeParam::idx_type,
                      typename TypeParam::idx_type>>
        e;

    t_tuple a{0, 0};
    t_tuple b{0, 1};
    t_tuple c{1, 0};
    t_tuple d{1, 1};
    e.push_back(t_tuple{1, 2});
    TypeParam tree(e, 4);

    k2_tree_test_nm::check_t_l(tree, {0, 1, 0, 0}, {0, 0, 1, 0});

    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {0, 1, 0, 0}, {0, 0, 1, 0});

    e.push_back(t_tuple{1, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {0, 1, 0, 0}, {0, 0, 1, 0});

    e.clear();
    e.push_back(t_tuple{0, 0});
    tree = TypeParam(e, 1);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0});

    e.push_back(t_tuple{0, 1});
    e.push_back(t_tuple{1, 0});
    e.push_back(t_tuple{1, 1});
    tree = TypeParam(e, 2);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 1, 1, 1});

    e.push_back(t_tuple{2, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {1, 0, 0, 1}, {1, 1, 1, 1, 1, 0, 0, 0});
}

TYPED_TEST(k2_tree_test_k_2, union_operation_test_happy_path)
{
    vector<vector<int>> mat({{1, 1, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 1, 1},
                             {0, 0, 1, 0}});
    // T 1 0 0 1
    // L 1 1 0 1  1 1 1 0
    TypeParam tree_A(mat);
    mat = vector<vector<int>>({{0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 1}});
    // T 0 0 0 1
    // L 0 0 0 1
    TypeParam tree_B(mat);

    TypeParam unionTree = tree_A.unionOp(tree_B);

    k2_tree_test_nm::check_t_l(unionTree, {1, 0, 0, 1}, {1, 1, 0, 1, 1, 1, 1, 1});
}

TYPED_TEST(k2_tree_test_k_2, empty_union_operation)
{
    vector<vector<int>> mat({{1, 1, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 1, 1},
                             {0, 0, 1, 0}});

    TypeParam tree_A(mat);
    TypeParam tree_B;

    TypeParam unionTree = tree_A.unionOp(tree_B);
    k2_tree_test_nm::check_t_l(unionTree, {1, 0, 0, 1}, {1, 1, 0, 1, 1, 1, 1, 0});

    unionTree = tree_B.unionOp(tree_A);

    k2_tree_test_nm::check_t_l(unionTree, {1, 0, 0, 1}, {1, 1, 0, 1, 1, 1, 1, 0});

    unionTree = tree_B.unionOp(tree_B);
    k2_tree_test_nm::assert_eq_tree(unionTree, tree_B);
}

TYPED_TEST(k2_tree_test_k_2, union_operation_test)
{
    vector<vector<int>> mat({{1, 1},
                             {0, 1}});
    TypeParam tree_A(mat);

    mat = vector<vector<int>>({{0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 1}});
    TypeParam tree_B(mat);

    ASSERT_THROW(tree_A.unionOp(tree_B), std::logic_error);
}

TYPED_TEST(k2_tree_test_k_2, edge_iterator_test)
{
    //forward iterator
    vector<vector<int>> mat({{1, 1, 0, 0},
                             {0, 0, 0, 0},
                             {0, 0, 1, 0},
                             {0, 0, 1, 0}});
    TypeParam tree(mat);
    auto edge_iterator = tree.edge_begin();
    ASSERT_EQ(std::get<0>(*edge_iterator), (size_t)0);
    ASSERT_EQ(std::get<1>(*edge_iterator), (size_t)0);

    // OPERATOR ASSIGNMENT
    auto another_edge_iterator = edge_iterator;
    ASSERT_EQ(std::get<0>(*another_edge_iterator), (size_t) 0);
    ASSERT_EQ(std::get<1>(*another_edge_iterator), (size_t) 0);

    // //OPERATOR INCREMENT
    // // ++edge_iterator; // also works
    edge_iterator++;
    ASSERT_EQ(std::get<0>(*edge_iterator), (size_t) 0);
    ASSERT_EQ(std::get<1>(*edge_iterator), (size_t) 1);
    edge_iterator++;
    ASSERT_EQ(std::get<0>(*edge_iterator), (size_t) 2);
    ASSERT_EQ(std::get<1>(*edge_iterator), (size_t) 2);
    edge_iterator++;
    ASSERT_EQ(std::get<0>(*edge_iterator), (size_t) 3);
    ASSERT_EQ(std::get<1>(*edge_iterator), (size_t) 2);
    edge_iterator++;

    //find last
    auto last = tree.edge_end();
    ASSERT_EQ(std::get<0>(*last), tree.get_number_nodes());
    ASSERT_EQ(std::get<1>(*last), tree.get_number_nodes());

    // OPERATOR EQUALS
    ASSERT_TRUE(edge_iterator == tree.edge_end());
    // OPERATOR INEQUALS
    ASSERT_TRUE(edge_iterator != tree.edge_begin());

    //Intensive test
    vector<vector<int>> mat2({{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 1, 0, 0, 1, 0, 0, 0, 0, 1},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}});
    TypeParam treeIntensive(mat2);

    auto it = treeIntensive.edge_begin();

    ASSERT_EQ(std::get<0>(*it), (size_t) 1);
    ASSERT_EQ(std::get<1>(*it), (size_t) 2);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 1);
    ASSERT_EQ(std::get<1>(*it), (size_t) 4);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 3);
    ASSERT_EQ(std::get<1>(*it), (size_t) 1);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 3);
    ASSERT_EQ(std::get<1>(*it), (size_t) 4);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 3);
    ASSERT_EQ(std::get<1>(*it), (size_t) 9);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 5);
    ASSERT_EQ(std::get<1>(*it), (size_t) 0);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 5);
    ASSERT_EQ(std::get<1>(*it), (size_t) 1);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 6);
    ASSERT_EQ(std::get<1>(*it), (size_t) 0);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 7);
    ASSERT_EQ(std::get<1>(*it), (size_t) 0);
    it++;
    ASSERT_EQ(std::get<0>(*it), (size_t) 9);
    ASSERT_EQ(std::get<1>(*it), (size_t) 0);

    //OPERATION SWAP
    swap(last, another_edge_iterator);
    ASSERT_EQ(std::get<0>(*last), (size_t) 0);
    ASSERT_EQ(std::get<1>(*last), (size_t) 0);
    ASSERT_EQ(std::get<0>(*another_edge_iterator), tree.get_number_nodes());
    ASSERT_EQ(std::get<1>(*another_edge_iterator), tree.get_number_nodes());
}

//     //OPERATOR DECREMENT
//     last = tree.edge_end();
//     last--;
//     ASSERT_EQ(std::get<0>(*last), (size_t) 3);
//     ASSERT_EQ(std::get<1>(*last), (size_t) 2);
//     last--;
//     ASSERT_EQ(std::get<0>(*last), (size_t) 2);
//     ASSERT_EQ(std::get<1>(*last), (size_t) 2);
//     last--;
//     ASSERT_EQ(std::get<0>(*last), (size_t) 0);
//     ASSERT_EQ(std::get<1>(*last), (size_t) 1);
//     last--; // does not go beyong the begin
//     ASSERT_EQ(std::get<0>(*last), (size_t) 0);
//     ASSERT_EQ(std::get<1>(*last), (size_t) 1);
// }

TYPED_TEST(k2_tree_test_k_2, node_iterator_test) {
    vector<vector<int>> mat({{0, 0, 0, 0},
                             {0, 0, 0, 0},
                             {0, 0, 1, 0},
                             {1, 1, 1, 0}});
    TypeParam tree(mat);
    auto node_iterator = tree.node_begin();
    ASSERT_EQ(*node_iterator, (size_t)0);

    node_iterator++;
    ASSERT_EQ(*node_iterator, (size_t)1);

    node_iterator++;
    ASSERT_EQ(*node_iterator, (size_t)2);

    node_iterator++;
    ASSERT_EQ(*node_iterator, (size_t)3);

    node_iterator++;
    ASSERT_EQ(*node_iterator, *tree.node_end());

    auto other_iterator = tree.node_begin();
    swap(other_iterator, node_iterator);
    ASSERT_EQ(*node_iterator, *tree.node_begin());
    ASSERT_EQ(*other_iterator, *tree.node_end());
}

TYPED_TEST(k2_tree_test_k_2, node_iterator_empty) {
    TypeParam empty_tree;
    ASSERT_EQ(empty_tree.node_begin(), empty_tree.node_end());
}

TYPED_TEST_CASE(k2_tree_test_k_3, k_3_implementations);

TYPED_TEST(k2_tree_test_k_3, build_from_matrix_test)
{
    vector<vector<int>> mat({{1, 1, 0, 0, 1},
                             {0, 1, 0, 0, 0},
                             {0, 0, 1, 1, 0},
                             {1, 1, 0, 1, 0},
                             {0, 0, 1, 0, 0}});

    TypeParam tree(mat);

    // cout << "ALGORIOTHM" << endl;
    // cout << all_of(tree.l().begin(), tree.l().end(), [](int i) { return i % 2 == 0; }) << endl;

    vector<unsigned> expected_t = {1, 1, 0, 1, 1, 0, 0, 0, 0};
    vector<unsigned> expected_l = {1, 1, 0, 0, 1, 0, 0, 0, 1,
                                   0, 1, 0, 0, 0, 0, 1, 0, 0,
                                   1, 1, 0, 0, 0, 1, 0, 0, 0,
                                   1, 0, 0, 0, 0, 0, 0, 0, 0};
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);

    mat = vector<vector<int>>({{1, 1, 1, 0},
                               {1, 0, 0, 0},
                               {0, 0, 0, 0},
                               {1, 1, 0, 0}});

    tree = TypeParam(mat);
    expected_t = {1, 0, 0, 1, 0, 0, 0, 0, 0};
    expected_l = {1, 1, 1, 1, 0, 0, 0, 0, 0,
                  1, 1, 0, 0, 0, 0, 0, 0, 0};
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);

    mat = vector<vector<int>>({{0, 0, 0},
                               {0, 0, 0},
                               {0, 0, 0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    // Size is minor than k:
    mat = vector<vector<int>>({{0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    mat = vector<vector<int>>({{1}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0, 0, 0, 0, 0, 0});

    mat = vector<vector<int>>({{1, 0},
                               {0, 1}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0, 1, 0, 0, 0, 0});

    // Size is a power of k:
    mat = vector<vector<int>>({{0, 0, 1, 0, 0, 0, 0, 0, 0},
                               {1, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 1, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 1, 0, 0, 0, 0, 0, 0}});
    tree = TypeParam(mat);
    expected_t = {1, 0, 0, 0, 0, 0, 1, 0, 0};
    expected_l = {0, 0, 1, 1, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1};
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);

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
    expected_t = {1, 1, 0, 1, 1, 0, 0, 0, 0,
                  1, 1, 0, 0, 0, 0, 0, 0, 1,
                  0, 0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 1, 0, 0, 0, 0, 0, 0,
                  1, 0, 0, 0, 0, 0, 0, 0, 0};

    expected_l = {0, 1, 0, 0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 1, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 0, 0,
                  1, 0, 1, 1, 0, 0, 0, 0, 0,
                  0, 1, 0, 1, 0, 0, 0, 0, 0};
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);
}

TYPED_TEST(k2_tree_test_k_3, build_from_edges_array)
{
    typedef std::tuple<typename TypeParam::idx_type,
                       typename TypeParam::idx_type>
        t_tuple;
    vector<std::tuple<typename TypeParam::idx_type,
                      typename TypeParam::idx_type>>
        e;

    e.push_back(t_tuple{1, 2});
    TypeParam tree(e, 4);

    k2_tree_test_nm::check_t_l(tree, {1, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 1, 0, 0, 0});

    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {}, {0, 0, 0, 0, 0, 1, 0, 0, 0});

    e.push_back(t_tuple{1, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {}, {0, 0, 0, 0, 0, 1, 0, 0, 0});

    e.clear();
    e.push_back(t_tuple{0, 0});
    tree = TypeParam(e, 1);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0, 0, 0, 0, 0, 0});

    e.push_back(t_tuple{0, 1});
    e.push_back(t_tuple{1, 0});
    e.push_back(t_tuple{1, 1});
    tree = TypeParam(e, 2);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 1, 0, 1, 1, 0, 0, 0, 0});

    e.clear();
    e.push_back(t_tuple{2, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {}, {0, 0, 0, 0, 0, 0, 0, 0, 1});
}

TYPED_TEST(k2_tree_test_k_3, union_operation_test)
{
    vector<vector<int>> mat({{1, 1, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 1, 1},
                             {0, 0, 1, 0}});

    TypeParam tree_A(mat);

    mat = vector<vector<int>>({{0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 0},
                               {0, 0, 0, 1}});
    TypeParam tree_B(mat);

    TypeParam unionTree = tree_A.unionOp(tree_B);
    k2_tree_test_nm::check_t_l(unionTree, {1, 1, 0, 1, 1, 0, 0, 0, 0}, {1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0});
}

TYPED_TEST(k2_tree_test_k_3, empty_union_operation)
{
    vector<vector<int>> mat({{1, 1, 0, 0, 1},
                             {0, 1, 0, 0, 0},
                             {0, 0, 1, 1, 0},
                             {1, 1, 0, 1, 0},
                             {0, 0, 1, 0, 0}});

    vector<unsigned> expected_t = {1, 1, 0, 1, 1, 0, 0, 0, 0};
    vector<unsigned> expected_l = {1, 1, 0, 0, 1, 0, 0, 0, 1,
                                   0, 1, 0, 0, 0, 0, 1, 0, 0,
                                   1, 1, 0, 0, 0, 1, 0, 0, 0,
                                   1, 0, 0, 0, 0, 0, 0, 0, 0};
    TypeParam tree_A(mat);
    TypeParam tree_B;

    TypeParam unionTree = tree_A.unionOp(tree_B);
    k2_tree_test_nm::check_t_l(unionTree, expected_t, expected_l);

    unionTree = tree_B.unionOp(tree_A);
    k2_tree_test_nm::check_t_l(unionTree, expected_t, expected_l);

    unionTree = tree_B.unionOp(tree_B);
    k2_tree_test_nm::assert_eq_tree(unionTree, tree_B);
}

TYPED_TEST_CASE(k2_tree_test, Implementations);

TYPED_TEST(k2_tree_test, edges_array_exhaustive)
{
    typedef std::tuple<typename TypeParam::idx_type,
                       typename TypeParam::idx_type>
        t_tuple;
    vector<std::tuple<typename TypeParam::idx_type,
                      typename TypeParam::idx_type>>
        e;
    e.push_back(t_tuple{5, 7});
    e.push_back(t_tuple{1, 2});
    e.push_back(t_tuple{3, 9});
    e.push_back(t_tuple{2, 2});
    e.push_back(t_tuple{3, 2});
    e.push_back(t_tuple{7, 5});
    e.push_back(t_tuple{1, 6});
    e.push_back(t_tuple{4, 8});
    e.push_back(t_tuple{4, 1});
    e.push_back(t_tuple{5, 2});

    TypeParam tree(e, 10);
    auto expected_neighbors = vector<vector<typename TypeParam::idx_type>>(10);
    expected_neighbors[0] = vector<typename TypeParam::idx_type>({});
    expected_neighbors[1] = vector<typename TypeParam::idx_type>({2, 6});
    expected_neighbors[2] = vector<typename TypeParam::idx_type>({2});
    expected_neighbors[3] = vector<typename TypeParam::idx_type>({2, 9});
    expected_neighbors[4] = vector<typename TypeParam::idx_type>({1, 8});
    expected_neighbors[5] = vector<typename TypeParam::idx_type>({2, 7});
    expected_neighbors[6] = vector<typename TypeParam::idx_type>({});
    expected_neighbors[7] = vector<typename TypeParam::idx_type>({5});
    expected_neighbors[8] = vector<typename TypeParam::idx_type>({});
    expected_neighbors[9] = vector<typename TypeParam::idx_type>({});
    for (unsigned i = 0; i < 10; i++)
    {
        auto actual_neighbors = tree.neigh(i);
        ASSERT_EQ(expected_neighbors[i].size(), actual_neighbors.size());
        for (unsigned j = 0; i < expected_neighbors[i].size(); i++)
            ASSERT_EQ(expected_neighbors[i][j], actual_neighbors[j]);
    }

    e.clear();
    e.push_back(t_tuple{0, 0});
    tree = TypeParam(e, 1);
    ASSERT_EQ(1u, tree.neigh(0).size());
    ASSERT_EQ(0u, tree.neigh(0)[0]);
}

TYPED_TEST(k2_tree_test, neighbors_test)
{
    vector<vector<int>> mat({{1, 1, 0, 0},
                             {0, 1, 0, 0},
                             {0, 0, 1, 1},
                             {0, 0, 1, 0}});

    TypeParam tree(mat);
    auto neigh_0 = tree.neigh(0);
    vector<unsigned> expected_neigh_0({0, 1});
    ASSERT_EQ(expected_neigh_0.size(), neigh_0.size());
    for (unsigned i = 0; i < neigh_0.size(); i++)
        ASSERT_EQ(expected_neigh_0[i], neigh_0[i]);

    auto neigh_3 = tree.neigh(3);
    vector<unsigned> expected_neigh_3({2});
    ASSERT_EQ(expected_neigh_3.size(), neigh_3.size());
    for (unsigned i = 0; i < neigh_3.size(); i++)
        ASSERT_EQ(expected_neigh_3[i], neigh_3[i]);

    mat = vector<vector<int>>({{1}});
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0[0]);
    ASSERT_EQ(1u, neigh_0.size());

    mat = vector<vector<int>>({{0, 0, 0},
                               {1, 0, 1},
                               {0, 1, 1}});
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0.size());

    auto neigh_1 = tree.neigh(1);
    auto expected_neigh_1 = vector<unsigned>({0, 2});
    ASSERT_EQ(expected_neigh_1.size(), neigh_1.size());
    for (unsigned i = 0; i < neigh_1.size(); i++)
        ASSERT_EQ(expected_neigh_1[i], neigh_1[i]);

    mat = vector<vector<int>>({{0, 0},
                               {0, 0}});
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0.size());
}

TYPED_TEST(k2_tree_test, reverse_neighbors_test)
{
    vector<vector<int>> mat({{1, 0, 0, 0, 1},
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

    for (unsigned i = 0; i < r_neigh_0.size(); i++)
        ASSERT_EQ(expected_r_neigh_0[i], r_neigh_0[i]);

    for (unsigned i = 0; i < r_neigh_2.size(); i++)
        ASSERT_EQ(expected_r_neigh_2[i], r_neigh_2[i]);

    mat = vector<vector<int>>({{0, 0},
                               {0, 0}});
    tree = TypeParam(mat);
    r_neigh_0 = tree.reverse_neigh(0);
    r_neigh_1 = tree.reverse_neigh(1);
    ASSERT_EQ(0u, r_neigh_0.size());
    ASSERT_EQ(0u, r_neigh_1.size());

    mat = vector<vector<int>>({{0, 1},
                               {1, 0}});
    tree = TypeParam(mat);
    r_neigh_0 = tree.reverse_neigh(0);
    expected_r_neigh_0 = vector<unsigned>({1});
    r_neigh_1 = tree.reverse_neigh(1);
    auto expected_r_neigh_1 = vector<unsigned>({0});

    ASSERT_EQ(expected_r_neigh_0.size(), r_neigh_0.size());
    ASSERT_EQ(expected_r_neigh_1.size(), r_neigh_1.size());
    for (unsigned i = 0; i < r_neigh_0.size(); i++)
        ASSERT_EQ(expected_r_neigh_0[i], r_neigh_0[i]);

    for (unsigned i = 0; i < r_neigh_1.size(); i++)
        ASSERT_EQ(expected_r_neigh_1[i], r_neigh_1[i]);
}

TYPED_TEST(k2_tree_test, range_test)
{
    static const vector<vector<int>> mat{
        //     0 1 2 3  4 5 6 7  8 9 10
        /* 0*/ {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        /* 1*/ {0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0},
        /* 2*/ {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        /* 3*/ {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

        /* 4*/ {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        /* 5*/ {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        /* 6*/ {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        /* 7*/ {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},

        /* 8*/ {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0},
        /* 9*/ {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1},
        /*10*/ {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0},
    };
    TypeParam tree(mat);

    auto result = tree.range(7, 9, 5, 9);
    static const decltype(result) expected{
        {7, 6},
        {8, 6},
        {8, 9},
        {9, 6},
        {9, 8}};
    sort(result.begin(), result.end());
    ASSERT_EQ(expected, result);

    auto all = tree.range(0, 10, 0, 10);
    static const decltype(all) expectedAll{
        {0, 1},
        {1, 2},
        {1, 3},
        {1, 4},
        {7, 6},
        {8, 6},
        {8, 9},
        {9, 6},
        {9, 8},
        {9, 10},
        {10, 6},
        {10, 9}};
    sort(all.begin(), all.end());
    ASSERT_EQ(expectedAll, all);

    auto single = tree.range(1, 1, 3, 3);
    static const decltype(single) expectedSingle{
        {1, 3}};
    ASSERT_EQ(expectedSingle, single);

    auto empty = tree.range(0, 6, 5, 10);
    static const decltype(empty) expectedEmpty{};
    ASSERT_EQ(expectedEmpty, empty);
}

TYPED_TEST(k2_tree_test, adj_test)
{
    vector<vector<int>> mat({{1, 0, 0, 0, 1},
                             {0, 0, 0, 0, 0},
                             {0, 0, 1, 1, 0},
                             {0, 0, 0, 0, 0},
                             {0, 0, 1, 0, 1}});

    auto tree = TypeParam(mat);
    ASSERT_TRUE(tree.adj(0, 0));
    ASSERT_TRUE(tree.adj(0, 4));
    ASSERT_FALSE(tree.adj(4, 0));
    ASSERT_TRUE(tree.adj(4, 4));
    ASSERT_FALSE(tree.adj(1, 1));
    ASSERT_TRUE(tree.adj(2, 2));
    ASSERT_TRUE(tree.adj(2, 3));

    mat = vector<vector<int>>({{0}});
    tree = TypeParam(mat);
    ASSERT_FALSE(tree.adj(0, 0));
    mat = vector<vector<int>>({{1}});
    tree = TypeParam(mat);
    ASSERT_TRUE(tree.adj(0, 0));
}

TYPED_TEST(k2_tree_test, serialize_test)
{
    vector<vector<int>> mat({{1, 0, 0, 0, 1},
                             {0, 0, 0, 0, 0},
                             {0, 0, 1, 1, 0},
                             {0, 0, 0, 0, 0},
                             {0, 0, 1, 0, 1}});

    auto tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);

    mat = vector<vector<int>>({{0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);

    tree = TypeParam();
    k2_tree_test_nm::check_serialize_load(tree);

    mat = vector<vector<int>>({{0, 0},
                               {0, 0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);

    mat = vector<vector<int>>({{1, 1},
                               {1, 1}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);
}

TYPED_TEST_CASE(k2_tree_test_marked, Implementations_bit_vector);
TYPED_TEST(k2_tree_test_marked, marked_edges)
{
    vector<vector<int>> mat({{1, 0, 0, 0, 1},
                             {0, 0, 0, 0, 0},
                             {0, 0, 1, 1, 0},
                             {0, 0, 0, 0, 0},
                             {0, 0, 1, 0, 1}});

    auto tree = TypeParam(mat);

    ASSERT_EQ(tree.get_marked_edges(), 0);
    ASSERT_EQ(tree.get_number_edges(), 6);
    ASSERT_TRUE(tree.erase(0, 0));
    ASSERT_EQ(tree.get_marked_edges(), 1);
    ASSERT_EQ(tree.get_number_edges(), 5);

    ASSERT_TRUE(tree.erase(0, 4));
    ASSERT_EQ(tree.get_marked_edges(), 2);
    ASSERT_EQ(tree.get_number_edges(), 4);

    ASSERT_FALSE(tree.erase(0, 4));
    ASSERT_EQ(tree.get_marked_edges(), 2);
    ASSERT_EQ(tree.get_number_edges(), 4);

    ASSERT_FALSE(tree.erase(1, 2));
    ASSERT_EQ(tree.get_marked_edges(), 2);
    ASSERT_EQ(tree.get_number_edges(), 4);
}

} // namespace

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
