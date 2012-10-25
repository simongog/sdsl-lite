#include "sdsl/wt_int.hpp"
#include "sdsl/util.hpp"
#include "sdsl/rrr_vector.hpp"
#include "sdsl/sd_vector.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;

template<class T>
class WtIntTest : public ::testing::Test
{
    protected:

        WtIntTest() {
            // You can do set-up work for each test here.
        }

        virtual ~WtIntTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            tmp_file = "tmp_wt_int_test_" + sdsl::util::to_string(sdsl::util::get_pid()) + "_";
			test_cases.push_back( sdsl::int_vector<>(1023,0,1) );
			test_cases.push_back( sdsl::int_vector<>(100023,0,1) );
			test_cases.push_back( sdsl::int_vector<>(64,0,2) );
			sdsl::int_vector<> iv(1000000,0,18);
    		sdsl::util::set_random_bits(iv, 17);
			test_cases.push_back( iv );
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }
        // Objects declared here can be used by all tests in the test case for Foo.
        std::string tmp_file;
        std::vector<sdsl::int_vector<> > test_cases;
};

using testing::Types;

typedef Types<
			  sdsl::wt_int<>,
     		  sdsl::wt_int<sdsl::int_vector<>, sdsl::rrr_vector<15> >,
     		  sdsl::wt_int<sdsl::int_vector<>, sdsl::rrr_vector<63> >
     		 > Implementations;

TYPED_TEST_CASE(WtIntTest, Implementations);

// TODO: test streaming operator

//! Test the parametrized constructor
TYPED_TEST(WtIntTest, Constructor)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		sdsl::int_vector<>& iv = this->test_cases[i];
		double iv_size = sdsl::util::get_size_in_mega_bytes(iv);
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
		sdsl::util::store_to_file(iv, tmp_file_name.c_str());
		{
			sdsl::int_vector_file_buffer<> buf(tmp_file_name.c_str());
			TypeParam wt(buf, buf.int_vector_size);
			std::cout << "compression = " << sdsl::util::get_size_in_mega_bytes(wt)/iv_size << std::endl;
			ASSERT_EQ(iv.size(), wt.size());
			for (size_type i=0; i < iv.size(); ++i) {
				ASSERT_EQ(iv[i], wt[i])<<i;
			}
		}
	}
}

TYPED_TEST(WtIntTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
        std::remove(tmp_file_name.c_str());
    }
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
