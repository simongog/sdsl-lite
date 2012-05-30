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
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }
        // Objects declared here can be used by all tests in the test case for Foo.
        std::string tmp_file;
};

using testing::Types;

typedef Types<
sdsl::wt_int<>,
     sdsl::wt_int<sdsl::int_vector<>, sdsl::rrr_vector<15> >,
     sdsl::wt_int<sdsl::int_vector<>, sdsl::rrr_vector<63> >
     > Implementations;

TYPED_TEST_CASE(WtIntTest, Implementations);

//! Test the parametrized constructor
TYPED_TEST(WtIntTest, Constructor)
{
    uint8_t width = 18;
    std::string suffix = "constructor";
    sdsl::int_vector<> iv(1000000,0,width);
    sdsl::util::set_random_bits(iv, 17);
    double iv_size = sdsl::util::get_size_in_mega_bytes(iv);

    sdsl::util::store_to_file(iv, (this->tmp_file+suffix).c_str());
    {
        sdsl::int_vector_file_buffer<> buf((this->tmp_file+suffix).c_str());
        TypeParam wt(buf);
        std::cout << "compression = " << sdsl::util::get_size_in_mega_bytes(wt)/iv_size << std::endl;
        ASSERT_EQ(iv.size(), wt.size());
        for (size_type i=0; i < iv.size(); ++i) {
            ASSERT_EQ(iv[i], wt[i])<<i;
        }
    }
    std::remove((this->tmp_file+suffix).c_str());
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
