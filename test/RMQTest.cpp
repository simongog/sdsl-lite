#include "sdsl/rmq_support.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <algorithm> // for std::min, random permutation
#include <stack>

namespace
{

typedef std::vector<uint64_t> tVUI;

ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;


void gen_increasing_vector(tVUI &v, tVUI::size_type len){
	v.resize(len);
	for( tVUI::size_type i=0; i < v.size(); ++i )
		v[i] = i;
}

template<class T>
class RMQTest : public ::testing::Test
{
    protected:

        RMQTest() {
            // You can do set-up work for each test here.
        }

        virtual ~RMQTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            tmp_file = "tmp_rmq_test_" + sdsl::util::to_string(sdsl::util::get_pid()) + "_";
			test_cases.push_back(tVUI(0)); // empty vector
			test_cases.push_back(tVUI(1,42)); // 1-element vector
			tVUI v;
			gen_increasing_vector(v, 1000000);
			test_cases.push_back(v); // increasing sequence
			srand(17);
			tVUI w;
			gen_increasing_vector(w, 1000000); 
			std::random_shuffle(w.begin(), w.end(), p_myrandom);
			test_cases.push_back(w); // random permutation
			std::random_shuffle(w.begin(), w.end(), p_myrandom);
			test_cases.push_back(w); // random permutation
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        template<class Rmq>
        std::string get_tmp_file_name(const Rmq& rmq, typename Rmq::size_type i) {
            std::locale loc;                 // the "C" locale
            const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc);
            std::string name = sdsl::util::demangle2(typeid(Rmq).name());
            uint64_t myhash = coll.hash(name.data(),name.data()+name.length());
            return tmp_file + sdsl::util::to_string(myhash) + "_" + sdsl::util::to_string(i);
        }

        template<class Rmq>
        bool load_rmq(Rmq& rmq, typename Rmq::size_type i) {
            return sdsl::util::load_from_file(rmq, get_tmp_file_name(rmq, i));
        }

        std::string tmp_file;
		std::vector<tVUI> test_cases;
};


using testing::Types;

typedef Types<  sdsl::rmq_succinct_sct<>,
				sdsl::rmq_succinct_sada<>
		     > Implementations;

TYPED_TEST_CASE(RMQTest, Implementations);


TYPED_TEST(RMQTest, ConstructAndStore)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		TypeParam rmq( &(this->test_cases[i]) );
        bool success = sdsl::util::store_to_file(rmq, this->get_tmp_file_name(rmq, i));
        ASSERT_EQ(success, true);
    }
}

// helper class for next test
class state{
	public:
	uint64_t l, r; // left and right border of interval
	uint64_t idx;  // index of the min value
	uint64_t min;  // min value in the interval
	state(uint64_t fl=0, uint64_t fr=0, uint64_t fidx = 0, uint64_t fmin=0) : 
		  l(fl), r(fr), idx(fidx), min(fmin) {}
};

//! Test range minimum queries
TYPED_TEST(RMQTest, RmqLoadAndQuery) {
	typedef typename TypeParam::size_type size_type;
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam rmq;
        ASSERT_EQ(this->load_rmq(rmq, i), true);
		ASSERT_EQ(rmq.size(), this->test_cases[i].size());
		if ( rmq.size() > 0 ) {
			std::stack<state> s;
			size_type idx = rmq(0, rmq.size()-1);
			ASSERT_TRUE( idx >= (size_type)0 ); ASSERT_TRUE( idx <= rmq.size() );
			s.push( state( 0, rmq.size()-1, idx, this->test_cases[i][idx] ) );
			while ( !s.empty() ){
				state st = s.top(); s.pop();
				if ( st.l < st.idx ){
					idx = rmq(st.l, st.idx-1);
					ASSERT_TRUE( idx >= st.l ); ASSERT_TRUE( idx <= st.idx-1 );
					ASSERT_TRUE( this->test_cases[i][idx] >= this->test_cases[i][st.idx] )
					<< "this->test_cases["<<i<<"]["<<idx<<"]="<< this->test_cases[i][idx]
					<< " < " << "this->test_cases["<<i<<"]["<<st.idx<<"]=" 
					<< this->test_cases[i][st.idx] << std::endl 
					<< "[" << st.l << "," << st.r << "]" << std::endl;
					s.push( state( st.l, st.idx-1, idx, this->test_cases[i][idx] ) );
				}
				if ( st.idx < st.r ){
					idx = rmq(st.idx+1, st.r);
					ASSERT_TRUE( idx >= st.idx+1 ); ASSERT_TRUE( idx <= st.r );
					ASSERT_TRUE( this->test_cases[i][idx] >= this->test_cases[i][st.idx] )
					<< "this->test_cases["<<i<<"]["<<idx<<"]="<< this->test_cases[i][idx]
					<< " < " << "this->test_cases["<<i<<"]["<<st.idx<<"]=" 
					<< this->test_cases[i][st.idx] << std::endl 
					<< "[" << st.l << "," << st.r << "]" << std::endl;
					s.push( state( st.idx+1, st.r, idx, this->test_cases[i][idx] ) );
				}
			}
		}
    }
}


TYPED_TEST(RMQTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam rmq;
        std::remove(this->get_tmp_file_name(rmq, i).c_str());
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
