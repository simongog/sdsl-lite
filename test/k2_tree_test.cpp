#include "sdsl/k2_tree.hpp"
#include "gtest/gtest.h"

#include <sstream>
#include <tuple>
#include <vector>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;

template<class T>
class k2_tree_test_k_2 : public ::testing::Test { };

template<class T>
class k2_tree_test_k_3 : public ::testing::Test { };

template<class T>
class k2_tree_test : public ::testing::Test { };

using testing::Types;

namespace k2_tree_test_nm
{
template<typename t_tree>
void check_t_l(t_tree& tree, vector<unsigned> expected_t,
               vector<unsigned> expected_l)
{
    ASSERT_EQ(expected_t.size(), tree.get_t().size());
    ASSERT_EQ(expected_l.size(), tree.get_l().size());
    for (unsigned i = 0; i < expected_t.size(); i++)
        ASSERT_EQ(expected_t[i], tree.get_t().get_int(i, 1));
    for (unsigned i = 0; i < expected_l.size(); i++)
        ASSERT_EQ(expected_l[i], tree.get_l().get_int(i, 1));
}

template<typename t_tree>
void check_serialize_load(t_tree& tree)
{
    auto unserialized_tree = t_tree();
    std::stringstream ss;
    tree.serialize(ss);
    unserialized_tree.load(ss);
    ASSERT_EQ(tree, unserialized_tree);
}
};

typedef Types<
k2_tree<2, bit_vector, rank_support_v<>>,
        k2_tree<2, bit_vector>
        > k_2_implementations;

typedef Types<
k2_tree<3, bit_vector, rank_support_v<>>,
        k2_tree<3, bit_vector>
        > k_3_implementations;

typedef Types<
k2_tree<2, bit_vector>,
        k2_tree<3, bit_vector>,
        k2_tree<7, bit_vector>,
        k2_tree<2, rrr_vector<63>>,
        k2_tree<3, rrr_vector<63>>,
        k2_tree<5, bit_vector, rank_support_v<>>,
        k2_tree<4, bit_vector, rank_support_v<>>
        > Implementations;

TYPED_TEST_CASE(k2_tree_test_k_2, k_2_implementations);

TYPED_TEST(k2_tree_test_k_2, build_from_matrix_test)
{
    vector<vector <int>> mat({{1, 1, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 1},
        {0, 0, 1, 0}
    });

    TypeParam tree(mat);
    vector<unsigned> expected_l = {1,1,0,1,1,1,1,0};
    k2_tree_test_nm::check_t_l(tree, {1, 0, 0 ,1}, expected_l);

    mat = vector<vector <int>> ({{0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
    });
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    mat = vector<vector<int>>({{0, 0},
        {0, 0}
    });
    tree = TypeParam(mat);
    ASSERT_TRUE(tree.get_t().empty());
    ASSERT_TRUE(tree.get_l().empty());

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
        {0, 1, 0}
    });
    tree = TypeParam(mat);
    expected_l = {0,0,0,1,1,0,0,0,0,1,0,0};
    k2_tree_test_nm::check_t_l(tree, {1, 1, 1 ,0}, expected_l);

    mat = vector<vector <int>>({{0, 0, 0},
        {1, 0, 1},
        {0, 1, 1}
    });
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
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0}
    });
    tree = TypeParam(mat);
    vector<unsigned> expected_t = {1, 0, 1, 1,  1, 1, 0, 1,  0, 1, 0, 0,
                                   1, 0, 0, 0,  1, 1, 0, 0,  1, 0, 0, 0,
                                   0, 0, 0, 1,  0, 1, 0, 1,  1, 1, 1, 0
                                  };

    expected_l = {0, 1, 0, 0,  0, 0, 1, 1,  0, 0, 1, 0,  0, 0, 1, 0,
                  1, 0, 1, 0,  1, 0, 0, 0,  0, 1, 1, 0,  0, 0, 1, 0,
                  0, 1, 0, 0
                 };
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);
}

TYPED_TEST(k2_tree_test_k_2, build_from_edges_array)
{
    typedef std::tuple<typename TypeParam::idx_type,
            typename TypeParam::idx_type> t_tuple;
    vector<std::tuple<typename TypeParam::idx_type,
           typename TypeParam::idx_type>> e;

    t_tuple a{0, 0};
    t_tuple b{0, 1};
    t_tuple c{1, 0};
    t_tuple d{1, 1};
    e.push_back(t_tuple {1, 2});
    TypeParam tree(e, 4);

    k2_tree_test_nm::check_t_l(tree, {0, 1, 0 ,0}, {0, 0, 1, 0});

    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {0, 1, 0 ,0}, {0, 0, 1, 0});

    e.push_back(t_tuple {1, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {0, 1, 0 ,0}, {0, 0, 1, 0});

    e.clear();
    e.push_back(t_tuple {0, 0});
    tree = TypeParam(e, 1);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0});

    e.push_back(t_tuple {0, 1});
    e.push_back(t_tuple {1, 0});
    e.push_back(t_tuple {1, 1});
    tree = TypeParam(e, 2);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 1, 1, 1});

    e.push_back(t_tuple {2, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {1, 0, 0, 1}, {1, 1, 1, 1,  1, 0, 0, 0});
}


TYPED_TEST_CASE(k2_tree_test_k_3, k_3_implementations);

TYPED_TEST(k2_tree_test_k_3, build_from_matrix_test)
{
    vector<vector <int>> mat({{1, 1, 0, 0, 1},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 1, 0},
        {1, 1, 0, 1, 0},
        {0, 0, 1, 0, 0}
    });

    TypeParam tree(mat);
    vector<unsigned> expected_t = {1, 1, 0, 1, 1, 0, 0, 0, 0};
    vector<unsigned> expected_l = {1, 1, 0, 0, 1, 0, 0, 0, 1,
                                   0, 1, 0, 0, 0, 0, 1, 0, 0,
                                   1, 1, 0, 0, 0, 1, 0, 0, 0,
                                   1, 0, 0, 0, 0, 0, 0, 0, 0
                                  };
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);

    mat = vector<vector<int>>({{1, 1, 1, 0},
        {1, 0, 0, 0},
        {0, 0, 0, 0},
        {1, 1, 0, 0}
    });

    tree = TypeParam(mat);
    expected_t = {1, 0, 0, 1, 0, 0, 0, 0, 0};
    expected_l = {1, 1, 1, 1, 0, 0, 0, 0, 0,
                  1, 1, 0, 0, 0, 0, 0, 0, 0
                 };
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);

    mat = vector<vector<int>>({{0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    });
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    // Size is minor than k:
    mat = vector<vector<int>>({{0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {});

    mat = vector<vector<int>>({{1}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0, 0, 0, 0, 0 ,0});

    mat = vector<vector<int>>({{1, 0},
        {0, 1}
    });
    tree = TypeParam(mat);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0, 1, 0, 0, 0 ,0});

    // Size is a power of k:
    mat = vector<vector<int>>({{0, 0, 1, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 0}
    });
    tree = TypeParam(mat);
    expected_t = {1, 0, 0, 0, 0, 0, 1, 0, 0};
    expected_l = {0, 0, 1, 1, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1
                 };
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
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0}
    });
    tree = TypeParam(mat);
    expected_t = {1, 1, 0, 1, 1, 0, 0, 0, 0,
                  1, 1, 0, 0, 0, 0, 0, 0, 1,
                  0, 0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 1, 0, 0, 0, 0, 0, 0,
                  1, 0, 0, 0, 0, 0, 0, 0, 0
                 };

    expected_l = {0, 1, 0, 0, 0, 1, 0, 0, 0,
                  0, 0, 0, 1, 1, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 0, 0,
                  1, 0, 1, 1, 0, 0, 0, 0, 0,
                  0, 1, 0, 1, 0, 0, 0, 0, 0
                 };
    k2_tree_test_nm::check_t_l(tree, expected_t, expected_l);
}


TYPED_TEST(k2_tree_test_k_3, build_from_edges_array)
{
    typedef std::tuple<typename TypeParam::idx_type,
            typename TypeParam::idx_type> t_tuple;
    vector<std::tuple<typename TypeParam::idx_type,
           typename TypeParam::idx_type>> e;

    e.push_back(t_tuple {1, 2});
    TypeParam tree(e, 4);

    k2_tree_test_nm::check_t_l(tree, {1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, 0, 0, 0});

    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {}, {0, 0, 0, 0, 0, 1, 0, 0, 0});

    e.push_back(t_tuple {1, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {}, {0, 0, 0, 0, 0, 1, 0, 0, 0});

    e.clear();
    e.push_back(t_tuple {0, 0});
    tree = TypeParam(e, 1);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 0, 0, 0, 0, 0, 0, 0, 0});

    e.push_back(t_tuple {0, 1});
    e.push_back(t_tuple {1, 0});
    e.push_back(t_tuple {1, 1});
    tree = TypeParam(e, 2);
    k2_tree_test_nm::check_t_l(tree, {}, {1, 1, 0, 1, 1, 0, 0, 0, 0});

    e.clear();
    e.push_back(t_tuple {2, 2});
    tree = TypeParam(e, 3);
    k2_tree_test_nm::check_t_l(tree, {}, {0, 0, 0, 0, 0, 0, 0, 0, 1});
}

TYPED_TEST_CASE(k2_tree_test, Implementations);

TYPED_TEST(k2_tree_test, edges_array_exhaustive)
{
    typedef std::tuple<typename TypeParam::idx_type,
            typename TypeParam::idx_type> t_tuple;
    vector<std::tuple<typename TypeParam::idx_type,
           typename TypeParam::idx_type>> e;
    e.push_back(t_tuple {5, 7});
    e.push_back(t_tuple {1, 2});
    e.push_back(t_tuple {3, 9});
    e.push_back(t_tuple {2, 2});
    e.push_back(t_tuple {3, 2});
    e.push_back(t_tuple {7, 5});
    e.push_back(t_tuple {1, 6});
    e.push_back(t_tuple {4, 8});
    e.push_back(t_tuple {4, 1});
    e.push_back(t_tuple {5, 2});

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
    for (unsigned i = 0; i < 10; i++) {
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
    vector<vector <int>> mat({{1, 1, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 1},
        {0, 0, 1, 0}
    });

    TypeParam tree(mat);
    auto neigh_0 = tree.neigh(0);
    vector<unsigned>expected_neigh_0({0, 1});
    ASSERT_EQ(expected_neigh_0.size(), neigh_0.size());
    for (unsigned i = 0; i < neigh_0.size(); i++)
        ASSERT_EQ(expected_neigh_0[i], neigh_0[i]);

    auto neigh_3 = tree.neigh(3);
    vector<unsigned>expected_neigh_3({2});
    ASSERT_EQ(expected_neigh_3.size(), neigh_3.size());
    for (unsigned i = 0; i < neigh_3.size(); i++)
        ASSERT_EQ(expected_neigh_3[i], neigh_3[i]);

    mat = vector<vector <int>>({{1}});
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0[0]);
    ASSERT_EQ(1u, neigh_0.size());

    mat = vector<vector <int>>({{0, 0, 0},
        {1, 0, 1},
        {0, 1, 1}
    });
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0.size());

    auto neigh_1 = tree.neigh(1);
    auto expected_neigh_1 = vector<unsigned>({0, 2});
    ASSERT_EQ(expected_neigh_1.size(), neigh_1.size());
    for (unsigned i = 0; i < neigh_1.size(); i++)
        ASSERT_EQ(expected_neigh_1[i], neigh_1[i]);

    mat = vector<vector <int>>({{0, 0},
        {0, 0}
    });
    tree = TypeParam(mat);
    neigh_0 = tree.neigh(0);
    ASSERT_EQ(0u, neigh_0.size());
}

TYPED_TEST(k2_tree_test, reverse_neighbors_test)
{
    vector<vector <int>> mat({{1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0},
        {0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1}
    });

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

    mat = vector<vector <int>>({{0, 0},
        {0, 0}
    });
    tree = TypeParam(mat);
    r_neigh_0 = tree.reverse_neigh(0);
    r_neigh_1 = tree.reverse_neigh(1);
    ASSERT_EQ(0u, r_neigh_0.size());
    ASSERT_EQ(0u, r_neigh_1.size());


    mat = vector<vector <int>>({{0, 1},
        {1, 0}
    });
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
		/* 0*/{0,1,0,0, 0,0,0,0, 0,0,0},
		/* 1*/{0,0,1,1, 1,0,0,0, 0,0,0},
		/* 2*/{0,0,0,0, 0,0,0,0, 0,0,0},
		/* 3*/{0,0,0,0, 0,0,0,0, 0,0,0},

		/* 4*/{0,0,0,0, 0,0,0,0, 0,0,0},
		/* 5*/{0,0,0,0, 0,0,0,0, 0,0,0},
		/* 6*/{0,0,0,0, 0,0,0,0, 0,0,0},
		/* 7*/{0,0,0,0, 0,0,1,0, 0,0,0},

		/* 8*/{0,0,0,0, 0,0,1,0, 0,1,0},
		/* 9*/{0,0,0,0, 0,0,1,0, 1,0,1},
		/*10*/{0,0,0,0, 0,0,1,0, 0,1,0},
	};
	TypeParam tree(mat);

	auto result = tree.range(7,9, 5,9);
	static const decltype(result) expected{
		{7,6},
		{8,6}, {8,9},
		{9,6}, {9,8}
	};
	sort(result.begin(), result.end());
	ASSERT_EQ(expected, result);

	auto all = tree.range(0,10, 0,10);
	static const decltype(all) expectedAll{
		{0,1},
		{1,2}, {1,3}, {1,4},
		{7,6},
		{8,6}, {8,9},
		{9,6}, {9,8}, {9,10},
		{10,6},{10,9}
	};
	sort(all.begin(), all.end());
	ASSERT_EQ(expectedAll, all);

	auto single = tree.range(1,1, 3,3);
	static const decltype(single) expectedSingle{
		{1,3}
	};
	ASSERT_EQ(expectedSingle, single);

	auto empty = tree.range(0,6, 5,10);
	static const decltype(empty) expectedEmpty{};
	ASSERT_EQ(expectedEmpty, empty);
}

TYPED_TEST(k2_tree_test, adj_test)
{
    vector<vector <int>> mat({{1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0},
        {0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1}
    });

    auto tree = TypeParam(mat);
    ASSERT_TRUE(tree.adj(0, 0));
    ASSERT_TRUE(tree.adj(0, 4));
    ASSERT_FALSE(tree.adj(4, 0));
    ASSERT_TRUE(tree.adj(4, 4));
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

TYPED_TEST(k2_tree_test, serialize_test)
{
    vector<vector <int>> mat({{1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0},
        {0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1}
    });

    auto tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);

    mat = vector<vector<int>>({{0}});
    tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);

    tree = TypeParam();
    k2_tree_test_nm::check_serialize_load(tree);

    mat = vector<vector<int>>({{0, 0},
        {0, 0}
    });
    tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);

    mat = vector<vector<int>>({{1, 1},
        {1, 1}
    });
    tree = TypeParam(mat);
    k2_tree_test_nm::check_serialize_load(tree);

}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
