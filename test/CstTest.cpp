#include "sdsl/suffixtrees.hpp"
#include "sdsl/test_index_performance.hpp"
#include "sdsl/util.hpp" // for store_to_file, load_to_file,...
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <locale>



namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::bit_vector bit_vector;

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
            test_cases.push_back("test_cases/crafted/100a.txt");
            test_cases.push_back("test_cases/small/faust.txt");
            test_cases.push_back("test_cases/small/zarathustra.txt");
            tmp_file = "tmp_cst_test_" + sdsl::util::to_string(sdsl::util::get_pid()) + "_";
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        std::vector<std::string> test_cases;
        std::string tmp_file;

        template<class Cst>
        std::string get_tmp_file_name(const Cst& cst) {
            std::locale loc;                 // the "C" locale
            const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc);
            std::string name = sdsl::util::demangle2(typeid(Cst).name());
            uint64_t myhash = coll.hash(name.data(),name.data()+name.length());
            return tmp_file + sdsl::util::to_string(myhash);
        }

        template<class Cst>
        bool load_cst(Cst& cst) {
            return sdsl::util::load_from_file(cst, get_tmp_file_name(cst).c_str());
        }
};

using testing::Types;

typedef Types<sdsl::cst_sct3<>,
        sdsl::cst_sada<>
        > Implementations;

TYPED_TEST_CASE(CstTest, Implementations);


TYPED_TEST(CstTest, CreateAndStoreTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        sdsl::util::verbose = false;
        construct_cst(this->test_cases[i], cst);
        bool success = sdsl::util::store_to_file(cst, this->get_tmp_file_name(cst).c_str());
        ASSERT_EQ(success, true);
    }
}

//! Test basic methods
TYPED_TEST(CstTest, BasicMethods)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(this->load_cst(cst), true);
        typedef typename TypeParam::node_type node_type;
        node_type r = cst.root(); // get root node
        // Size of the subtree rooted at r should the size of the suffix array
        ASSERT_EQ(cst.leaves_in_the_subtree(r), cst.csa.size());
        // Check leaf methods
        for (size_type i=0; i < cst.csa.size(); ++i) {
            ASSERT_EQ(cst.is_leaf(cst.ith_leaf(i+1)), true);
        }
    }
}


//! Test the id method
/*
TYPED_TEST(CstTest, IdMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
		ASSERT_EQ( this->load_cst(cst), true);
		// doing a depth first traversal through the tree to count the nodes
		typedef typename TypeParam::const_iterator const_iterator;
		typedef typename TypeParam::node_type node_type;
		size_type node_count=0;
		for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it){
			if ( it.visit() == 1 ){
				++node_count;
			}
		}
		// counted nodes should be equal to nodes
		ASSERT_EQ( node_count, cst.nodes() );
		// check if the id method is working
		bit_vector marked(cst.nodes(), 0);
		for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it){
			if ( it.visit() == 1 ){
				++node_count;
				node_type v = *it;
				size_type id = cst.id( v );
				ASSERT_EQ( marked[id], 0 );
				marked[id] = 1;
				ASSERT_EQ( cst.inv_id( cst.id(v) ), v );
			}
		}

    }
}
*/


template<class Cst>
typename Cst::node_type naive_lca(const Cst& cst, typename Cst::node_type v, typename Cst::node_type w)
{
    size_type steps = 0;
    while (v != w  and steps < cst.csa.size()) {
        if (cst.depth(v) > cst.depth(w)) {
            v = cst.parent(v);
        } else {
            w = cst.parent(w);
        }
        steps++;
    }
    return v;
}


TYPED_TEST(CstTest, LcaMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(this->load_cst(cst), true);
        uint64_t mask;
        uint8_t log_m = 14;
        // create m/2 pairs of positions in [0..cst.csa.size()-1]
        typedef typename TypeParam::node_type node_type;
        sdsl::int_vector<64> rnd_pos = sdsl::get_rnd_positions(log_m, mask, cst.csa.size());
        // test for random sampled nodes
        for (size_type i=0; i < rnd_pos.size()/2; ++i) {
            // get two children
            node_type v = cst.ith_leaf(rnd_pos[2*i]+1);
            node_type w = cst.ith_leaf(rnd_pos[2*i+1]+1);
            // calculate lca
            node_type z = naive_lca(cst, v, w);
            ASSERT_EQ(cst.lca(v, w), z);
        }
        // test for regular sampled nodes
        for (size_type i=cst.csa.size()/2, g=100; i+g < cst.csa.size()/2+100*g; ++i) {
            // get two children
            node_type v = cst.ith_leaf(i+1);
            node_type w = cst.ith_leaf(i+g+1);
            // calculate lca
            node_type z = naive_lca(cst, v, w);
            ASSERT_EQ(cst.lca(v, w), z);
        }
    }
}

//! Test the node method
TYPED_TEST(CstTest, NodeMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(this->load_cst(cst), true);
        // doing a depth first traversal through the tree to count the nodes
        typedef typename TypeParam::const_iterator const_iterator;
        typedef typename TypeParam::node_type node_type;
        for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it) {
            if (it.visit() == 1) {
                node_type v = *it;
                size_type d = cst.depth(v);
                size_type lb = cst.lb(v), rb = cst.rb(v);
                ASSERT_EQ(cst.node(lb, rb, d), v);
            }
        }
    }
}


//! Test the bottom-up iterator
TYPED_TEST(CstTest, BottomUpIterator)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(this->load_cst(cst), true);
        // doing a bottom-up traversal of the tree
        // TODO: implement
    }
}

TYPED_TEST(CstTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        std::remove(this->get_tmp_file_name(cst).c_str());
    }
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

