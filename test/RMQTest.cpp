#include "sdsl/rmq_support.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <algorithm> // for std::min

namespace
{

typedef std::vector<uint64_t> tVUI;


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
			tVUI v;
			gen_increasing_vector(v, 1000000);
			test_cases.push_back(v);
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
            return sdsl::util::load_from_file(rmq, get_tmp_file_name(rmq, i).c_str());
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
        bool success = sdsl::util::store_to_file(rmq, this->get_tmp_file_name(rmq, i).c_str());
        ASSERT_EQ(success, true);
    }
}

//! Test range minimum queries
TYPED_TEST(RMQTest, RmqLoadAndQuery) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(this->load_rmq(csa, i), true);
    }
}


TYPED_TEST(RMQTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        std::remove(this->get_tmp_file_name(csa, i).c_str());
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
