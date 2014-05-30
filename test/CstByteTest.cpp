#include "CstHelper.hpp"
#include "sdsl/suffix_trees.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <random>

using namespace sdsl;
using namespace std;

namespace
{

typedef int_vector<>::size_type size_type;
typedef bit_vector bit_vector;

tMSS  test_case_file_map;
string test_file;
string temp_file;
string temp_dir;
bool in_memory;

template<class T>
class CstByteTest : public ::testing::Test { };

using testing::Types;

typedef Types<
cst_sct3<>,
         cst_sada<>,
         cst_sct3<cst_sct3<>::csa_type, lcp_bitcompressed<>>,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_tree2<>>,
         cst_sada<cst_sada<>::csa_type, lcp_dac<>>,
         cst_sada<cst_sada<>::csa_type, lcp_vlc<>>,
         cst_sada<cst_sada<>::csa_type, lcp_byte<>>,
         cst_sada<cst_sada<>::csa_type, lcp_support_tree2<>, bp_support_gg<>>,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_tree<>, bp_support_gg<>>,
         cst_sada<cst_sada<>::csa_type, lcp_support_tree<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_sada<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_wt<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_tree<>, bp_support_g<> >,
         cst_sct3<csa_bitcompressed<>, lcp_bitcompressed<> >
         > Implementations;

TYPED_TEST_CASE(CstByteTest, Implementations);


TYPED_TEST(CstByteTest, CreateAndStoreTest)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    TypeParam cst;
    ASSERT_TRUE(cst.empty());
    cache_config config(false, temp_dir, util::basename(test_file));
    construct(cst, test_file, config, 1);
    test_case_file_map = config.file_map;
    ASSERT_TRUE(store_to_file(cst, temp_file));
    TypeParam cst2;
    cst2 = cst;
    ASSERT_EQ(cst.size(), cst2.size());
    ASSERT_TRUE(cst.size() <= TypeParam::max_size());
}

//! Test the swap method
TYPED_TEST(CstByteTest, SwapMethod)
{
    TypeParam cst1;
    ASSERT_TRUE(load_from_file(cst1, temp_file));
    size_type n = cst1.size();
    TypeParam cst2;
    ASSERT_EQ((size_type)0, cst2.size());
    cst1.swap(cst2);
    ASSERT_EQ((size_type)0, cst1.size());
    ASSERT_EQ(n, cst2.size());
    ASSERT_EQ(n, cst2.csa.size());
    bit_vector mark((size_type)0, cst2.size());
    check_node_method(cst2);
}

//! Test the move method
TYPED_TEST(CstByteTest, MoveMethod)
{
    TypeParam cst1;
    ASSERT_TRUE(load_from_file(cst1, temp_file));
    size_type n = cst1.size();
    TypeParam cst2 = std::move(cst1);
    ASSERT_EQ(n, cst2.size());
    ASSERT_EQ(n, cst2.csa.size());
    bit_vector mark((size_type)0, cst2.size());
    check_node_method(cst2);
}


//! Test the node method
TYPED_TEST(CstByteTest, NodeMethod)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    // doing a depth first traversal through the tree to count the nodes
    check_node_method(cst);
}

//! Test basic methods
TYPED_TEST(CstByteTest, BasicMethods)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    typedef typename TypeParam::node_type node_type;
    node_type r = cst.root(); // get root node
    // Size of the subtree rooted at r should the size of the suffix array
    ASSERT_EQ(cst.csa.size(), cst.size(r));
    // Check leaf methods
    for (size_type i=0; i < cst.csa.size(); ++i) {
        ASSERT_TRUE(cst.is_leaf(cst.select_leaf(i+1)));
    }
}

//! Test suffix array access
TYPED_TEST(CstByteTest, SaAccess)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    sdsl::int_vector<> sa;
    sdsl::load_from_file(sa, test_case_file_map[sdsl::conf::KEY_SA]);
    size_type n = sa.size();
    ASSERT_EQ(n, cst.csa.size());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(sa[j], cst.csa[j])<<" j="<<j;
    }
}

//! Test suffix array access after move
TYPED_TEST(CstByteTest, MoveSaAccess)
{
    TypeParam cst_load;
    ASSERT_TRUE(load_from_file(cst_load, temp_file));
    TypeParam cst = std::move(cst_load);
    sdsl::int_vector<> sa;
    sdsl::load_from_file(sa, test_case_file_map[sdsl::conf::KEY_SA]);
    size_type n = sa.size();
    ASSERT_EQ(n, cst.csa.size());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(sa[j], cst.csa[j])<<" j="<<j;
    }
}

//! Test BWT access
TYPED_TEST(CstByteTest, BwtAccess)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    sdsl::int_vector<8> bwt;
    sdsl::load_from_file(bwt, test_case_file_map[sdsl::conf::KEY_BWT]);
    size_type n = bwt.size();
    ASSERT_EQ(n, cst.csa.bwt.size());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(bwt[j], cst.csa.bwt[j])<<" j="<<j;
    }
}

//! Test LCP access
TYPED_TEST(CstByteTest, LcpAccess)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    sdsl::int_vector<> lcp;
    sdsl::load_from_file(lcp, test_case_file_map[sdsl::conf::KEY_LCP]);
    size_type n = lcp.size();
    ASSERT_EQ(n, cst.lcp.size());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(lcp[j], cst.lcp[j])<<" j="<<j;
    }
}

//! Test LCP access after move
TYPED_TEST(CstByteTest, MoveLcpAccess)
{
    TypeParam cst_load;
    ASSERT_TRUE(load_from_file(cst_load, temp_file));
    TypeParam cst = std::move(cst_load);
    sdsl::int_vector<> lcp;
    sdsl::load_from_file(lcp, test_case_file_map[sdsl::conf::KEY_LCP]);
    size_type n = lcp.size();
    ASSERT_EQ(n, cst.lcp.size());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(lcp[j], cst.lcp[j])<<" j="<<j;
    }
}

//! Test the id and inverse id method
TYPED_TEST(CstByteTest, IdMethod)
{
    TypeParam cst;
    // test empty iterator
    ASSERT_EQ(cst.begin(), cst.end());

    ASSERT_TRUE(load_from_file(cst, temp_file));
    // doing a depth first traversal through the tree to count the nodes
    typedef typename TypeParam::const_iterator const_iterator;
    typedef typename TypeParam::node_type node_type;
    size_type node_count=0;
    for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it) {
        if (it.visit() == 1) {
            ++node_count;
        }
    }
    // counted nodes should be equal to nodes
    ASSERT_EQ(node_count, cst.nodes());
    // check if the id method is working
    bit_vector marked(cst.nodes(), 0);
    for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it) {
        if (it.visit() == 1) {
            ++node_count;
            node_type v = *it;
            size_type id = cst.id(v);
            ASSERT_EQ(0, marked[id]);
            marked[id] = 1;
            ASSERT_EQ(v, cst.inv_id(cst.id(v)));
        }
    }
}

TYPED_TEST(CstByteTest, SelectChild)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    if (cst.size() > 1) {
        ASSERT_EQ(cst.csa.sigma, cst.degree(cst.root()));
        size_type lb = 0;
        for (size_type i=1; i <= cst.csa.sigma; ++i) {
            auto v = cst.select_child(cst.root(), i);
            ASSERT_EQ(lb, cst.lb(v));
            lb = cst.rb(v)+1;
        }
        ASSERT_EQ(cst.rb(cst.root()), lb-1);

        size_type i=1;
        for (auto v  : cst.children(cst.root())) {
            ASSERT_TRUE(i <= cst.degree(cst.root()));
            ASSERT_EQ(cst.select_child(cst.root(),i), v) << i << "!";
            ++i;
        }
    } else if (cst.size() == 1) {
        ASSERT_EQ(1U, cst.csa.sigma);
        ASSERT_EQ(0U, cst.degree(cst.root()));
    }
}

TYPED_TEST(CstByteTest, MoveSelectChild)
{
    TypeParam cst_load;
    ASSERT_TRUE(load_from_file(cst_load, temp_file));
    TypeParam cst = std::move(cst_load);
    if (cst.size() > 1) {
        ASSERT_EQ(cst.csa.sigma, cst.degree(cst.root()));
        size_type lb = 0;
        for (size_type i=1; i <= cst.csa.sigma; ++i) {
            auto v = cst.select_child(cst.root(), i);
            ASSERT_EQ(lb, cst.lb(v));
            lb = cst.rb(v)+1;
        }
        ASSERT_EQ(cst.rb(cst.root()), lb-1);

        size_type i=1;
        for (auto v  : cst.children(cst.root())) {
            ASSERT_TRUE(i <= cst.degree(cst.root()));
            ASSERT_EQ(cst.select_child(cst.root(),i), v) << i << "!";
            ++i;
        }
    } else if (cst.size() == 1) {
        ASSERT_EQ(1U, cst.csa.sigma);
        ASSERT_EQ(0U, cst.degree(cst.root()));
    }
}



TYPED_TEST(CstByteTest, SelectLeafAndSn)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    for (size_type i=0; i < std::min(cst.csa.size(), (size_type)100); ++i) {
        ASSERT_EQ(cst.csa[i], cst.sn(cst.select_leaf(i+1)));
    }
}


TYPED_TEST(CstByteTest, NodeDepth)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    auto v = cst.root();
    ASSERT_EQ((size_type)0, cst.node_depth(v));
    for (size_type i=1; i<=10 and !cst.is_leaf(v); ++i) {
        v = cst.select_child(v, 2);
        ASSERT_EQ(i, cst.node_depth(v));
    }
}


TYPED_TEST(CstByteTest, Child)
{
    TypeParam cst;
    typedef typename TypeParam::char_type char_type;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    if (cst.size() > 1) {
        std::set<char_type> char_set;
        ASSERT_EQ(cst.csa.sigma, cst.degree(cst.root()));
        for (size_type i=0; i < cst.csa.sigma; ++i) {
            auto c = cst.csa.comp2char[i];
            char_set.insert(c);
            auto v = cst.select_child(cst.root(), i+1);
            auto w = cst.child(cst.root(), c);
            ASSERT_EQ(v, w);
            if (cst.is_leaf(v)) {
                ASSERT_EQ(cst.root(), cst.select_child(v, c));
            }
        }
        for (size_type i=0; i < 256; ++i) {
            char_type c = (char_type)i;
            if (char_set.find(c) == char_set.end()) {
                ASSERT_EQ(cst.root(), cst.child(cst.root(), c));
            }
        }
    }
}

TYPED_TEST(CstByteTest, Edge)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));

    int_vector<8> data;
    ASSERT_TRUE(load_vector_from_file(data, test_file, 1));

    if (cst.csa.size() > 0) {
        auto v = cst.select_leaf(cst.csa.isa[0]+1);
        size_type max_depth = std::min(cst.depth(v), (size_type)20);
        for (size_type i=0; i<max_depth; ++i) {
            ASSERT_EQ(data[i], cst.edge(v, i+1))<<" i="<<i<<" v="<<v;
        }
        v = cst.parent(v);
        max_depth = std::min(max_depth, cst.depth(v));
        for (size_type i=0; i<max_depth; ++i) {
            ASSERT_EQ(data[i], cst.edge(v, i+1))<<" i="<<i<<" v="<<v;
        }
    }
}

TYPED_TEST(CstByteTest, LeftmostRightmostLeaf)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    if (cst.size() > 0) {
        auto v = cst.select_leaf(cst.size()/2+1);
        while (true) {
            auto v_l = cst.leftmost_leaf(v);
            auto v_r = cst.rightmost_leaf(v);
            ASSERT_TRUE(cst.is_leaf(v_l));
            ASSERT_TRUE(cst.is_leaf(v_r));
            ASSERT_EQ(cst.lb(v), cst.lb(v_l));
            ASSERT_EQ(cst.rb(v), cst.rb(v_r));
            if (v == cst.root())
                break;
            v = cst.parent(v);
        }
    }
}

TYPED_TEST(CstByteTest, SuffixAndWeinerLink)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    ASSERT_EQ(cst.root(),cst.sl(cst.root()));

    if (cst.size() > 0) {
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, cst.size()-1);
        auto dice = bind(distribution, rng);

        for (size_type i=0; i<100; ++i) {
            auto v = cst.select_leaf(dice()+1);
            auto c = cst.edge(v, 1);
            ASSERT_EQ(v, cst.wl(cst.sl(v), c));
            for (size_type j=0; j<5; ++j) {
                v = cst.parent(v);
                if (cst.root() == v)
                    break;
                c = cst.edge(v, 1);
                ASSERT_EQ(v, cst.wl(cst.sl(v), c));
            }
        }
    }
}

TYPED_TEST(CstByteTest, LcaMethod)
{
    TypeParam cst;
    ASSERT_TRUE(load_from_file(cst, temp_file));
    uint64_t mask;
    uint8_t log_m = 6;
    // create m/2 pairs of positions in [0..cst.csa.size()-1]
    typedef typename TypeParam::node_type node_type;
    int_vector<64> rnd_pos = util::rnd_positions<int_vector<64>>(log_m, mask, cst.csa.size());
    // test for random sampled nodes
    for (size_type i=0; i < rnd_pos.size()/2; ++i) {
        // get two children
        node_type v = cst.select_leaf(rnd_pos[2*i]+1);
        node_type w = cst.select_leaf(rnd_pos[2*i+1]+1);
        // calculate lca
        node_type z = naive_lca(cst, v, w);
        ASSERT_EQ(z, cst.lca(v, w));
    }
    // test for regular sampled nodes
    size_type g = std::max(cst.csa.size()/30, (size_type)5);
    for (size_type i=cst.csa.size()/2; i+g < cst.csa.size(); ++i) {
        // get two children
        node_type v = cst.select_leaf(i+1);
        node_type w = cst.select_leaf(i+g+1);
        // calculate lca
        node_type z = naive_lca(cst, v, w);
        node_type u = cst.lca(v, w);
        ASSERT_EQ(z, u) << " naive_lca is "
                        << naive_lca(cst, v, w, true) << endl;
    }
}

//! Test the bottom-up iterator
TYPED_TEST(CstByteTest, BottomUpIterator)
{
//    TypeParam cst;
//    ASSERT_TRUE(load_from_file(cst, temp_file));
//    doing a bottom-up traversal of the tree
//    TODO: implement
}

TYPED_TEST(CstByteTest, DeleteTest)
{
    sdsl::remove(temp_file);
    util::delete_all_files(test_case_file_map);
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 4) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " test_file temp_file tmp_dir [in-memory]" << endl;
        cout << " (1) Generates a CST out of test_file; stores it in temp_file." << endl;
        cout << "     Temporary files (SA/BWT/LCP/TEXT) are stored in tmp_dir." << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    temp_file = argv[2];
    temp_dir  = argv[3];
    in_memory    = argc > 4;
    if (in_memory) {
        temp_dir = "@";
        int_vector<8> data;
        load_vector_from_file(data, test_file, 1);
        test_file = ram_file_name(test_file);
        store_to_plain_array<uint8_t>(data, test_file);
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}

