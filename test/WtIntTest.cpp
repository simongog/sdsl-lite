#include "sdsl/wt_int.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;

// The fixture for testing class int_vector.
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
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }
        // Objects declared here can be used by all tests in the test case for Foo.
        std::string tmp_file;
};


//! Test the parametrized constructor
TEST_F(WtIntTest, Constructor)
{
	uint8_t width = 18;
	std::string suffix = "constructor";
    sdsl::int_vector<> iv(10000000,0,width);
	sdsl::util::set_random_bits(iv, 17);

	sdsl::util::store_to_file(iv, (tmp_file+suffix).c_str() );
	{
		sdsl::int_vector_file_buffer<> buf( (tmp_file+suffix).c_str() ); 
		sdsl::wt_int<> wt(buf);
		ASSERT_EQ( iv.size(), wt.size() );
		for (size_type i=0; i < iv.size(); ++i) {
			ASSERT_EQ(iv[i], wt[i]);
		}
	}
	std::remove( (tmp_file+suffix).c_str() );
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
