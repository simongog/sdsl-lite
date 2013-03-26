#include "sdsl/int_vector.hpp"
#include "sdsl/bits.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;

// The fixture for testing class int_vector.
class BitMagicTest : public ::testing::Test
{
    protected:

        BitMagicTest() {
            // You can do set-up work for each test here.
        }

        virtual ~BitMagicTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test case for Foo.
};

//! Test the default constructor
TEST_F(BitMagicTest, Rank)
{
    for (size_type i=0; i<64; ++i) {
        ASSERT_EQ(sdsl::bits::cnt(1ULL<<i), (uint64_t)1);
    }
}

//! Test the parametrized constructor
TEST_F(BitMagicTest, Select)
{
    for (size_type i=0; i<64; ++i) {
        ASSERT_EQ(sdsl::bits::sel(1ULL<<i, 1), i);
        ASSERT_EQ(sdsl::bits::selv(1ULL<<i, 1), i);
        ASSERT_EQ(sdsl::bits::sel3(1ULL<<i, 1), i);
    }
    for (size_type i=0; i < 10000000; ++i) {
        uint64_t x = (((size_type)rand())<<32) + rand();
        for (size_type j=0, ones=0; j<64; ++j) {
            if ((x >> j)&1) {
                ++ones;
                ASSERT_EQ(sdsl::bits::sel(x, ones), j);
                ASSERT_EQ(sdsl::bits::selv(x, ones), j);
                ASSERT_EQ(sdsl::bits::sel3(x, ones), j);
            }
        }
    }
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
