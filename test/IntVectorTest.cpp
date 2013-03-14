#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::int_vector<>::value_type value_type;

// The fixture for testing class int_vector.
class IntVectorTest : public ::testing::Test
{
    protected:

        IntVectorTest() {}

        virtual ~IntVectorTest() {}

        virtual void SetUp() {
            vec_sizes.push_back(0);
            vec_sizes.push_back(64);
            vec_sizes.push_back(65);
            vec_sizes.push_back(127);
            vec_sizes.push_back(128);
            for (size_type i=0; i<128; ++i) {
                vec_sizes.push_back(rand() % 100000);
            }
            for (size_type i=0; i < 10; ++i) {
                vec_sizes.push_back(rand()%10000000);
            }
        }

        virtual void TearDown() {}

        std::vector<size_type> vec_sizes; // different sizes for the vectors
};

//! Test the default constructor
TEST_F(IntVectorTest, DefaultConstruct)
{
    {
        sdsl::int_vector<> iv;
        ASSERT_EQ((size_type)0, iv.size());			   // size should be 0
        ASSERT_EQ((uint8_t)64, iv.width());   // default width of each element should be 64 bits
    }
    {
        sdsl::bit_vector bv;
        ASSERT_EQ((size_type)0, bv.size());
        ASSERT_EQ(1, bv.width());  // for a bit vector the width should be 1 bits
    }
    {
        sdsl::int_vector<8> iv;
        ASSERT_EQ((size_type)0, iv.size());
        ASSERT_EQ(8, iv.width());   // default width of each element should be 8 bits
    }
    {
        sdsl::int_vector<16> iv;
        ASSERT_EQ((size_type)0, iv.size());
        ASSERT_EQ(16, iv.width());   // default width of each element should be 16 bits
    }
    {
        sdsl::int_vector<32> iv;
        ASSERT_EQ((size_type)0, iv.size());
        ASSERT_EQ(32, iv.width());   // default width of each element should be 32 bits
    }
    {
        sdsl::int_vector<64> iv;
        ASSERT_EQ((size_type)0, iv.size());
        ASSERT_EQ(64, iv.width());   // default width of each element should be 64 bits
    }
}

//! Test the parametrized constructor
TEST_F(IntVectorTest, Constructor)
{
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        sdsl::int_vector<> iv(vec_sizes[i]);
        ASSERT_EQ(vec_sizes[i], iv.size());
        for (size_type j=0; j < iv.size(); ++j) { // should be initialized with 0s
            ASSERT_EQ((size_type)0, iv[j]);
        }
    }
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        const size_type val = rand();
        sdsl::int_vector<> iv(vec_sizes[i], val);
        ASSERT_EQ(vec_sizes[i], iv.size());
        for (size_type j=0; j < iv.size(); ++j) { // each element should be initialized with val
            ASSERT_EQ(val, iv[j]);
        }
    }
}

TEST_F(IntVectorTest, Swap)
{
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        const size_type val = rand();
        sdsl::int_vector<> iv(vec_sizes[i], val);
        {
            sdsl::int_vector<> tmp;
            ASSERT_EQ((size_type)0, tmp.size());
            tmp.swap(iv);
            ASSERT_EQ((size_type)0, iv.size());
            ASSERT_EQ(vec_sizes[i], tmp.size());
            for (size_type j=0; j < tmp.size(); ++j) {
                ASSERT_EQ(val, tmp[j]);
            }
        }
    }
}

TEST_F(IntVectorTest, AssignElement)
{
    for (unsigned char w=1; w <= 64; ++w) { // for each possible width
        sdsl::int_vector<> iv(100000, 0, w);
        for (size_type i=0; i < iv.size(); ++i) {
            value_type x = rand() & sdsl::bit_magic::Li1Mask[w];
            iv[i] = x;
            ASSERT_EQ(x, iv[i]);
        }
    }
}

TEST_F(IntVectorTest, AddAndSub)
{
    for (unsigned char w=1; w <= 64; ++w) { // for each possible width
        sdsl::int_vector<> iv(100000, 0, w);
        sdsl::util::set_random_bits(iv);
        for (size_type i=0; i < iv.size(); ++i) {
            value_type x = iv[i];
            value_type y = rand() & sdsl::bit_magic::Li1Mask[w];
            iv[i] += y;
            ASSERT_EQ((x+y)&sdsl::bit_magic::Li1Mask[w], iv[i]);
            iv[i] -= y;
            ASSERT_EQ( x & sdsl::bit_magic::Li1Mask[w], iv[i]);
            iv[i]++;
            ASSERT_EQ((x+1)&sdsl::bit_magic::Li1Mask[w], iv[i]);
            iv[i]--;
            ASSERT_EQ( x & sdsl::bit_magic::Li1Mask[w], iv[i]);
            ++iv[i];
            ASSERT_EQ((x+1)&sdsl::bit_magic::Li1Mask[w], iv[i]);
            --iv[i];
            ASSERT_EQ( x & sdsl::bit_magic::Li1Mask[w], iv[i]);
        }
    }
}

TEST_F(IntVectorTest, SerializeAndLoad)
{
    sdsl::int_vector<> iv(1000000);
    for (size_type i=0; i<iv.size(); ++i)
        iv[i] = rand();
    std::string file_name = "/tmp/int_vector";
    sdsl::util::store_to_file(iv, file_name);
    sdsl::int_vector<> iv2;
    sdsl::util::load_from_file(iv2, file_name);
    ASSERT_EQ(iv.size(), iv2.size());
    for (size_type i=0; i<iv.size(); ++i)
        ASSERT_EQ(iv[i], iv2[i]);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
