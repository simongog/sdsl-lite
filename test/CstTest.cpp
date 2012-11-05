#include "sdsl/suffixtrees.hpp"
#include "sdsl/lcp.hpp"
#include "sdsl/test_index_performance.hpp"
#include "sdsl/util.hpp" // for store_to_file, load_to_file,...
#include "sdsl/config.hpp" // for CMAKE_SOURCE_DIR
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <sstream>

using namespace sdsl;

namespace
{

typedef int_vector<>::size_type size_type;
typedef bit_vector bit_vector;

template<class T>
class CstTest : public ::testing::Test
{
    protected:

        CstTest() {
            // You can do set-up work for each test here
        }

        virtual ~CstTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            string test_cases_dir = string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/test_cases";
            test_cases.push_back(test_cases_dir + "/crafted/100a.txt");
            test_cases.push_back(test_cases_dir + "/small/faust.txt");
            test_cases.push_back(test_cases_dir + "/small/zarathustra.txt");
            tmp_file = "tmp_cst_test_" + util::to_string(util::get_pid()) + "_";
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        std::vector<std::string> test_cases;
        std::string tmp_file;

        template<class Cst>
        std::string get_tmp_file_name(const Cst& cst, size_type i) {
            return tmp_file + util::class_to_hash(cst) + "_" + util::basename(test_cases[i].c_str());
        }

        template<class Cst>
        bool load_cst(Cst& cst, size_type i) {
            return util::load_from_file(cst, get_tmp_file_name(cst, i).c_str());
        }
};

using testing::Types;

typedef Types<
		cst_sct3<cst_sada<>::csa_type, lcp_bitcompressed<> >,
         cst_sada<cst_sada<>::csa_type, lcp_dac<> >,
         cst_sada<cst_sada<>::csa_type, lcp_vlc<> >,
         cst_sada<cst_sada<>::csa_type, lcp_support_tree2<>, bp_support_gg<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_tree<>, bp_support_gg<> >,
         cst_sct3<>,
         cst_sada<>,
         cst_sada<cst_sada<>::csa_type, lcp_support_tree<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_tree2<> >,
         cst_sada<cst_sada<>::csa_type, lcp_dac<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_sada<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_tree<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_wt<> >,
         cst_sada<cst_sada<>::csa_type, lcp_support_sada<> >,
         cst_sada<cst_sada<>::csa_type, lcp_support_tree<> >,
         cst_sada<cst_sada<>::csa_type, lcp_support_tree2<> >,
         cst_sada<cst_sada<>::csa_type, lcp_wt<> >,
         cst_sct3<cst_sct3<>::csa_type, lcp_support_tree<>, bp_support_g<> >,
		 cst_sct3<csa_bitcompressed<>, lcp_bitcompressed<> >, 
         cst_sada<cst_sada<>::csa_type, lcp_dac<>, bp_support_g<> >
         > Implementations;

TYPED_TEST_CASE(CstTest, Implementations);


TYPED_TEST(CstTest, CreateAndStoreTest)
{
    util::verbose=true;
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        util::verbose = false;
        construct_cst(this->test_cases[i], cst);
        bool success = util::store_to_file(cst, this->get_tmp_file_name(cst, i).c_str());
        ASSERT_EQ(true, success);
    }
}


template<class tCst>
void check_node_method(const tCst& cst)
{
    typedef typename tCst::const_iterator const_iterator;
    typedef typename tCst::node_type node_type;
    for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it) {
        if (it.visit() == 1) {
            node_type v = *it;
            size_type d = cst.depth(v);
            size_type lb = cst.lb(v), rb = cst.rb(v);
            ASSERT_EQ(v, cst.node(lb, rb, d));
        }
    }
}



//! Test the swap method
TYPED_TEST(CstTest, SwapMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst1;
        ASSERT_EQ(true, this->load_cst(cst1, i));
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
}


//! Test the node method
TYPED_TEST(CstTest, NodeMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(true, this->load_cst(cst, i));
        // doing a depth first traversal through the tree to count the nodes
        check_node_method(cst);
    }
}

//! Test basic methods
TYPED_TEST(CstTest, BasicMethods)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(true, this->load_cst(cst, i));
        typedef typename TypeParam::node_type node_type;
        node_type r = cst.root(); // get root node
        // Size of the subtree rooted at r should the size of the suffix array
        ASSERT_EQ(cst.csa.size(), cst.leaves_in_the_subtree(r));
        // Check leaf methods
        for (size_type i=0; i < cst.csa.size(); ++i) {
            ASSERT_EQ(true, cst.is_leaf(cst.select_leaf(i+1)));
        }
    }
}


//! Test the id and inverse id method
TYPED_TEST(CstTest, IdMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(true, this->load_cst(cst, i));
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
}



template<class Cst>
std::string format_node(const Cst& cst, const typename Cst::node_type& v)
{
    std::stringstream ss;
    ss << cst.depth(v) << "-["<<cst.lb(v)<<","<<cst.rb(v)<<"]";
    return ss.str();
}

template<class Cst>
typename Cst::node_type naive_lca(const Cst& cst, typename Cst::node_type v, typename Cst::node_type w, bool output=false)
{
    size_type steps = 0;
    while (v != w  and steps < cst.csa.size()) {
        if (cst.depth(v) > cst.depth(w)) {
            v = cst.parent(v);
            if (output) {
                std::cout << "v="<<format_node(cst, v) << std::endl;
            }
        } else {
            w = cst.parent(w);
            if (output) {
                std::cout << "w="<<format_node(cst, v) << std::endl;
            }
        }
        steps++;
    }
    return v;
}


TYPED_TEST(CstTest, LcaMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(true, this->load_cst(cst, i));
        uint64_t mask;
        uint8_t log_m = 14;
        // create m/2 pairs of positions in [0..cst.csa.size()-1]
        typedef typename TypeParam::node_type node_type;
        int_vector<64> rnd_pos = get_rnd_positions(log_m, mask, cst.csa.size());
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
        for (size_type i=cst.csa.size()/2, g=100; i+g < std::min(cst.csa.size(), cst.csa.size()/2+100*g); ++i) {
            // get two children
            node_type v = cst.select_leaf(i+1);
            node_type w = cst.select_leaf(i+g+1);
            // calculate lca
            node_type z = naive_lca(cst, v, w);
            node_type u = cst.lca(v, w);
            if (u != z) {
                std::cout << "v="<<v<<" w="<<w<< std::endl;
                std::cout << "u="<<u<<" z="<<z<<std::endl;
                std::cout << "--------------"<<std::endl;
                std::cout << "v="<<format_node(cst, v)<<" w="<<format_node(cst, w)<< std::endl;
                std::cout << "u="<<format_node(cst, u)<<" z="<<format_node(cst, z)<<std::endl;
                naive_lca(cst, v, w, true);
            }
            ASSERT_EQ(z, u);
        }
    }
}





//! Test the bottom-up iterator
TYPED_TEST(CstTest, BottomUpIterator)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
//        TypeParam cst;
//        ASSERT_EQ(this->load_cst(cst, i), true);
        // doing a bottom-up traversal of the tree
        // TODO: implement
    }
}

TYPED_TEST(CstTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        std::remove(this->get_tmp_file_name(cst, i).c_str());
    }
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

