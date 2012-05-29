#include "sdsl/suffixtrees.hpp"
#include "sdsl/lcp.hpp"
#include "sdsl/test_index_performance.hpp"
#include "sdsl/util.hpp" // for store_to_file, load_to_file,...
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <locale>
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
            test_cases.push_back("test_cases/crafted/100a.txt");
            test_cases.push_back("test_cases/small/faust.txt");
            test_cases.push_back("test_cases/small/zarathustra.txt");
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
            std::locale loc;                 // the "C" locale
            const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc);
            std::string name = util::demangle2(typeid(Cst).name());
            uint64_t myhash = coll.hash(name.data(),name.data()+name.length());
            return tmp_file + util::to_string(myhash) + "_" + util::basename(test_cases[i].c_str());
        }

        template<class Cst>
        bool load_cst(Cst& cst, size_type i) {
            return util::load_from_file(cst, get_tmp_file_name(cst, i).c_str());
        }
};

using testing::Types;

typedef Types<cst_sct3<>,
        cst_sada<>,
        cst_sct3<cst_sct3<>::csa_type, lcp_support_sada<> >,
        cst_sct3<cst_sct3<>::csa_type, lcp_support_tree<> >,
        cst_sct3<cst_sct3<>::csa_type, lcp_support_tree2<> >,
        cst_sct3<cst_sct3<>::csa_type, lcp_wt<> >
        > Implementations;

TYPED_TEST_CASE(CstTest, Implementations);


TYPED_TEST(CstTest, CreateAndStoreTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        util::verbose = false;
        construct_cst(this->test_cases[i], cst);
        bool success = util::store_to_file(cst, this->get_tmp_file_name(cst, i).c_str());
        ASSERT_EQ(success, true);
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
            ASSERT_EQ(cst.node(lb, rb, d), v);
        }
    }
}

//! Test the swap method
TYPED_TEST(CstTest, SwapMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst1;
        ASSERT_EQ(this->load_cst(cst1, i), true);
        size_type n = cst1.size();
        TypeParam cst2;
        ASSERT_EQ(cst2.size(), 0);
        cst1.swap(cst2);
        ASSERT_EQ(cst1.size(), 0);
        ASSERT_EQ(cst2.size(), n);
        ASSERT_EQ(cst2.csa.size(), n);
        bit_vector mark(cst2.size(), 0);
        /*		for(size_type i=0; i<cst2.size(); ++i){
        			size_type x = cst2.csa[i];
        			ASSERT_EQ( mark[x], 0 );
        			mark[x] = 1;
        		}
        */
        check_node_method(cst2);
    }
}

//! Test basic methods
TYPED_TEST(CstTest, BasicMethods)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(this->load_cst(cst, i), true);
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
        ASSERT_EQ(this->load_cst(cst, i), true);
        uint64_t mask;
        uint8_t log_m = 14;
        // create m/2 pairs of positions in [0..cst.csa.size()-1]
        typedef typename TypeParam::node_type node_type;
        int_vector<64> rnd_pos = get_rnd_positions(log_m, mask, cst.csa.size());
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
        for (size_type i=cst.csa.size()/2, g=100; i+g < std::min(cst.csa.size(), cst.csa.size()/2+100*g); ++i) {
            // get two children
            node_type v = cst.ith_leaf(i+1);
            node_type w = cst.ith_leaf(i+g+1);
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
            ASSERT_EQ(u, z);
        }
    }
}


//! Test the node method
TYPED_TEST(CstTest, NodeMethod)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam cst;
        ASSERT_EQ(this->load_cst(cst, i), true);
        // doing a depth first traversal through the tree to count the nodes
        check_node_method(cst);
        /*		typedef typename TypeParam::const_iterator const_iterator;
                typedef typename TypeParam::node_type node_type;
                for (const_iterator it = cst.begin(), end = cst.end(); it != end; ++it) {
                    if (it.visit() == 1) {
                        node_type v = *it;
                        size_type d = cst.depth(v);
                        size_type lb = cst.lb(v), rb = cst.rb(v);
                        ASSERT_EQ(cst.node(lb, rb, d), v);
                    }
                }
        */
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

