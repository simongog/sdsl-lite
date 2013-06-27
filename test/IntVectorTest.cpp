#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <random>
#include <algorithm>

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
            rng.seed((uint32_t)time(NULL));
            vec_sizes.push_back(64);
            vec_sizes.push_back(0);
            vec_sizes.push_back(65);
            vec_sizes.push_back(127);
            vec_sizes.push_back(128);
            for (size_type i=0; i<128; ++i) {
                vec_sizes.push_back(rng() % 100000);
            }
            for (size_type i=0; i < 10; ++i) {
                vec_sizes.push_back(rng()%10000000);
            }
        }

        virtual void TearDown() {}

        std::mt19937_64 rng; // random number generator
        std::vector<size_type> vec_sizes; // different sizes for the vectors
};

template<class t_iv>
void test_Constructors(bool special, uint8_t template_width, size_type constructor_size, uint8_t constructor_width)
{
    std::mt19937_64 rng((uint32_t)time(NULL));
    {
        // Constructor without argument
        t_iv iv;
        ASSERT_EQ((size_type)0, iv.size()); // default size should be 0
        ASSERT_EQ((uint8_t)template_width, iv.width()); // verify default width of each element
    }
    {
        // Constructor with one argument
        t_iv iv(constructor_size);
        ASSERT_EQ(constructor_size, iv.size());
        ASSERT_EQ(template_width, iv.width());
        for (size_type j=0; j < iv.size(); ++j) { // should be initialized with 0s
            ASSERT_EQ((typename t_iv::value_type)0, (typename t_iv::value_type)iv[j]);
        }
    }
    {
        // Constructor with two arguments
        size_type expected_val = rng();
        t_iv iv(constructor_size, expected_val);
        ASSERT_EQ(constructor_size, iv.size());
        ASSERT_EQ(template_width, iv.width());
        expected_val &= sdsl::bits::lo_set[iv.width()];
        for (size_type j=0; j < iv.size(); ++j) { // should be initialized with expected_val
            ASSERT_EQ(expected_val, (size_type)iv[j]);
        }
    }
    {
        // Constructor with three arguments
        size_type expected_val = rng();
        t_iv iv(constructor_size, expected_val, constructor_width);
        ASSERT_EQ(constructor_size, iv.size());
        if (special) {
            ASSERT_EQ(template_width, iv.width());
        } else {
            ASSERT_EQ(constructor_width, iv.width());
        }
        expected_val &= sdsl::bits::lo_set[iv.width()];
        for (size_type j=0; j < iv.size(); ++j) { // should be initialized with expected_val
            ASSERT_EQ(expected_val, (size_type)iv[j]);
        }
    }
}

//! Test Constructors
TEST_F(IntVectorTest, Constructors)
{
    std::vector<size_type> vec_widths = {15, 32, 55}; // different widths for the vectors
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        for (size_type j=0; j<vec_widths.size(); ++j) {
            // unspecialized
            test_Constructors<sdsl::int_vector<>   >(false, 64, vec_sizes[i], vec_widths[j]);
            test_Constructors<sdsl::int_vector<3>  >(false,  3, vec_sizes[i], vec_widths[j]);
            test_Constructors<sdsl::int_vector<31> >(false, 31, vec_sizes[i], vec_widths[j]);
            // specialized
            test_Constructors<sdsl::bit_vector     >(true,  1, vec_sizes[i], vec_widths[j]);
            test_Constructors<sdsl::int_vector<8>  >(true,  8, vec_sizes[i], vec_widths[j]);
            test_Constructors<sdsl::int_vector<16> >(true, 16, vec_sizes[i], vec_widths[j]);
            test_Constructors<sdsl::int_vector<32> >(true, 32, vec_sizes[i], vec_widths[j]);
            test_Constructors<sdsl::int_vector<64> >(true, 64, vec_sizes[i], vec_widths[j]);
        }
    }
}

TEST_F(IntVectorTest, Swap)
{
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        const size_type val = rng();
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
            value_type x = rng() & sdsl::bits::lo_set[w];
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
            value_type y = rng() & sdsl::bits::lo_set[w];
            iv[i] += y;
            ASSERT_EQ((x+y)&sdsl::bits::lo_set[w], iv[i]);
            iv[i] -= y;
            ASSERT_EQ(x & sdsl::bits::lo_set[w], iv[i]);
            iv[i]++;
            ASSERT_EQ((x+1)&sdsl::bits::lo_set[w], iv[i]);
            iv[i]--;
            ASSERT_EQ(x & sdsl::bits::lo_set[w], iv[i]);
            ++iv[i];
            ASSERT_EQ((x+1)&sdsl::bits::lo_set[w], iv[i]);
            --iv[i];
            ASSERT_EQ(x & sdsl::bits::lo_set[w], iv[i]);
        }
    }
}

TEST_F(IntVectorTest, STL)
{
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        sdsl::int_vector<> iv(vec_sizes[i]);
        ASSERT_EQ(vec_sizes[i], iv.size());
        auto cnt = iv.size();
for (auto x : iv) {
            x = --cnt;
        }
        std::sort(iv.begin(), iv.end());
        sdsl::int_vector<>::value_type last = 0;
for (const auto& x : iv) {
            ASSERT_TRUE(x >= last);
            last = x;
        }
    }
}

TEST_F(IntVectorTest, SerializeAndLoad)
{
    sdsl::int_vector<> iv(1000000);
    for (size_type i=0; i<iv.size(); ++i)
        iv[i] = rng();
    std::string file_name = "/tmp/int_vector";
    sdsl::store_to_file(iv, file_name);
    sdsl::int_vector<> iv2;
    sdsl::load_from_file(iv2, file_name);
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
