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
            std::mt19937_64 rng;
            {
                std::uniform_int_distribution<uint64_t> distribution(0, 100000);
                auto dice = bind(distribution, rng);
                for (size_type i=0; i<128; ++i) {
                    vec_sizes.push_back(dice());
                }
            }
            {
                std::uniform_int_distribution<uint64_t> distribution(0, 10000000);
                auto dice = bind(distribution, rng);
                for (size_type i=0; i < 10; ++i) {
                    vec_sizes.push_back(dice());
                }
            }
        }

        virtual void TearDown() {}

        std::vector<size_type> vec_sizes = {0,64,65,127,128}; // different sizes for the vectors
};

template<class t_iv>
void test_Constructors(uint8_t template_width, size_type constructor_size, uint8_t constructor_width)
{
    std::mt19937_64 rng;
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
        if (iv.fixed_int_width == 0) {
            ASSERT_EQ(constructor_width, iv.width());
        } else {
            ASSERT_EQ(template_width, iv.width());
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
    for (auto size : vec_sizes) {
        if (size<1000) {                                // Test only for short sizes,
            for (uint8_t width=1; width<=64; ++width) { // but for all possible widths
                // unspecialized
                test_Constructors<sdsl::int_vector<>   >(64, size, width);
                test_Constructors<sdsl::int_vector<3>  >(3, size, width);
                test_Constructors<sdsl::int_vector<31> >(31, size, width);
                // specialized
                test_Constructors<sdsl::bit_vector     >(1, size, width);
                test_Constructors<sdsl::int_vector<8>  >(8, size, width);
                test_Constructors<sdsl::int_vector<16> >(16, size, width);
                test_Constructors<sdsl::int_vector<32> >(32, size, width);
                test_Constructors<sdsl::int_vector<64> >(64, size, width);
            }
        }
    }
}

TEST_F(IntVectorTest, Swap)
{
    std::mt19937_64 rng;
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

template<class t_iv>
void test_AssignAndModifyElement(uint64_t size, uint8_t width)
{
    std::mt19937_64 rng;
    t_iv iv(size, 0, width);
    for (size_type i=0; i<iv.size(); ++i) {
        value_type exp_v = rng(), tmp = rng();

        iv[i] = exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]);
        iv[i] += tmp;
        exp_v += tmp;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]);
        iv[i] -= tmp;
        exp_v -= tmp;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]++);
        exp_v++;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]--);
        exp_v--;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]);
        ++exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], ++iv[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]);
        --exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], --iv[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[width], iv[i]);
    }
}

template<>
void test_AssignAndModifyElement<sdsl::bit_vector>(uint64_t size, uint8_t width)
{
    std::mt19937_64 rng(13);
    std::uniform_int_distribution<uint64_t> distribution(0, 9);

    sdsl::bit_vector bv(size, 0, width);
    for (size_type i=0; i<bv.size(); ++i) {
        value_type exp_v = distribution(rng);
        bv[i] = exp_v;
        ASSERT_EQ((bool)exp_v, bv[i]);
    }
    bv.flip();

    rng.seed(13); // To get the same values
    for (size_type i=0; i<bv.size(); ++i) {
        value_type exp_v = !distribution(rng);
        ASSERT_EQ((bool)exp_v, bv[i]);
    }
}

TEST_F(IntVectorTest, AssignAndModifyElement)
{
    // unspecialized vector for each possible width
    for (uint8_t width=1; width <= 64; ++width) {
        test_AssignAndModifyElement< sdsl::int_vector<> >(100000, width);
    }
    // specialized vectors
    test_AssignAndModifyElement<sdsl::bit_vector     >(100000,  1);
    test_AssignAndModifyElement<sdsl::int_vector< 8> >(100000,  8);
    test_AssignAndModifyElement<sdsl::int_vector<16> >(100000, 16);
    test_AssignAndModifyElement<sdsl::int_vector<32> >(100000, 32);
    test_AssignAndModifyElement<sdsl::int_vector<64> >(100000, 64);
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
    std::mt19937_64 rng;
    sdsl::int_vector<> iv(1000000);
    for (size_type i=0; i<iv.size(); ++i)
        iv[i] = rng();
    std::string file_name = "tmp/int_vector";
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
