#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector_buffer.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::int_vector<>::value_type value_type;

// The fixture for testing class int_vector.
class IntVectorBufferTest : public ::testing::Test
{
    protected:

        IntVectorBufferTest() {}

        virtual ~IntVectorBufferTest() {}

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

        std::vector<size_type> vec_sizes = {0, 64, 65, 127, 128}; // different sizes for the vectors
};



template<class t_T>
void test_constructors(size_type template_width, size_type constructor_width, size_type exp_width)
{
    std::string file_name = "tmp/int_vector_buffer";
    {
        // Default constructor
        t_T ivb;
        ASSERT_EQ("", ivb.filename());                    // filename should be empty
        ASSERT_EQ((size_type)0, ivb.size());                         // size should be 0
        ASSERT_EQ((size_type)0, ivb.buffersize());                   // buffersize should be 0
        ASSERT_EQ((uint8_t)template_width, ivb.width());  // default width of each element should be template_width bits
    }
    {
        // Constructor with 2 Parameters
        t_T ivb(file_name, false);
        ASSERT_EQ(file_name, ivb.filename());                               // filename should be file_name
        ASSERT_EQ((size_type)0, ivb.size());                                           // size should be 0
        ASSERT_EQ((uint8_t)template_width, ivb.width());                    // default width of each element should be template_width bits
        ASSERT_LT((ivb.buffersize()-(1024*1024)), ivb.width());             // actual buffersize() is not more than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
    }
    {
        // Constructor with 4 Parameters
        size_type buffersize = 1024*1024;
        t_T ivb(file_name, false, buffersize, constructor_width);
        ASSERT_EQ(file_name, ivb.filename());                   // filename should be file_name
        ASSERT_EQ((size_type)0, ivb.size());                               // size should be 0
        ASSERT_EQ((uint8_t)exp_width, ivb.width());             // default width of each element should be i bits
        ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
    }
}

//! Test constructors
TEST_F(IntVectorBufferTest, Constructors)
{
    for (size_type width=1; width<=64; ++width) {
        test_constructors< sdsl::int_vector_buffer<> >(64, width, width);
        test_constructors< sdsl::int_vector_buffer<1> >(1, width, 1);
        test_constructors< sdsl::int_vector_buffer<8> >(8, width, 8);
        test_constructors< sdsl::int_vector_buffer<13> >(13, width, 13);
        test_constructors< sdsl::int_vector_buffer<16> >(16, width, 16);
        test_constructors< sdsl::int_vector_buffer<32> >(32, width, 32);
        test_constructors< sdsl::int_vector_buffer<64> >(64, width, 64);
    }
}


template<class t_T>
void test_buffersize(std::vector<size_type>& vec_sizes)
{
    std::string file_name = "tmp/int_vector_buffer";
    // Test constructor with different buffersizes
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        t_T ivb(file_name, false, vec_sizes[i]);
        ASSERT_EQ(file_name, ivb.filename());                     // filename should be file_name
        ASSERT_EQ((size_type)0, ivb.size());                                 // size should be 0
        ASSERT_LT((ivb.buffersize()-vec_sizes[i]), ivb.width());  // actual buffersize() is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
    }
    {
        // Test buffersize() method
        t_T ivb(file_name, false);
        for (size_type i=0; i < vec_sizes.size(); ++i) {
            ivb.buffersize(vec_sizes[i]);
            ASSERT_LT((ivb.buffersize()-vec_sizes[i]), ivb.width());  // actual buffersize() is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
    }
}

//! Test the buffersize in constructor and buffersize()
TEST_F(IntVectorBufferTest, Buffersize)
{
    test_buffersize< sdsl::int_vector_buffer<> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<1> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<8> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<13> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<16> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<32> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<64> >(vec_sizes);
}


template<class t_T>
void test_assign_and_modify(size_type width=1)
{
    std::mt19937_64 rng(13), rng2;
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;
    t_T ivb(file_name, false, buffersize, width);
    for (size_type i=0; i < 100000; ++i) {
        value_type exp_v = rng();
        ivb[i] = exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
    }
    rng.seed(13); // To get the same values
    for (size_type i=0; i < 100000; ++i) {
        value_type exp_v = rng(), tmp = rng2();
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ivb[i] += tmp;
        exp_v += tmp;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ivb[i] -= tmp;
        exp_v -= tmp;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]++);
        exp_v++;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]--);
        exp_v--;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ++exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)++ivb[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        --exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)--ivb[i]);
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
    }
}

template<>
void test_assign_and_modify<sdsl::int_vector_buffer<1>>(size_type width)
{
    std::mt19937_64 rng(13);
    std::uniform_int_distribution<uint64_t> distribution(0, 9);

    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;
    sdsl::int_vector_buffer<1> ivb(file_name, false, buffersize, width);
    for (size_type i=0; i < 100000; ++i) {
        value_type exp_v = distribution(rng);
        ivb[i] = exp_v;
        ASSERT_EQ((bool)exp_v, (size_type)ivb[i]);
    }
    rng.seed(13); // To get the same values
    for (size_type i=0; i < 100000; ++i) {
        value_type exp_v = distribution(rng);
        ASSERT_EQ((bool)exp_v, (size_type)ivb[i]);
    }
}

TEST_F(IntVectorBufferTest, AssignAndModifyElement)
{
    for (size_type width=1; width<=64; ++width) { // for each possible width
        test_assign_and_modify< sdsl::int_vector_buffer<> >(width);
    }
    test_assign_and_modify< sdsl::int_vector_buffer<1> >();
    test_assign_and_modify< sdsl::int_vector_buffer<8> >();
    test_assign_and_modify< sdsl::int_vector_buffer<13> >();
    test_assign_and_modify< sdsl::int_vector_buffer<16> >();
    test_assign_and_modify< sdsl::int_vector_buffer<32> >();
    test_assign_and_modify< sdsl::int_vector_buffer<64> >();
}

template<class t_T>
void compare(size_type width=1)
{
    std::mt19937_64 rng(13);
    std::string file_name_1 = "tmp/int_vector_buffer_1";
    std::string file_name_2 = "tmp/int_vector_buffer_2";
    size_type buffersize = 1024;
    t_T ivb1(file_name_1, false, buffersize, width);
    t_T ivb2(file_name_2, false, buffersize, width);
    for (size_type i=0; i < 1000; ++i) {
        ivb1[i] = rng() & sdsl::bits::lo_set[ivb1.width()];
        ivb2[i] = ivb1[i];
        ASSERT_EQ(true, ivb1[i]==ivb2[i]);
        ASSERT_EQ(false, ivb1[i]!=ivb2[i]);
        ASSERT_EQ((size_type)ivb1[i], (size_type)ivb2[i]);
    }
    for (size_type i=0; i < 1000; ++i) {
        ivb1[i] = rng() & sdsl::bits::lo_set[ivb1.width()];
        ivb2[i] = rng() & sdsl::bits::lo_set[ivb1.width()];
        size_type v1 = ivb1[i];
        size_type v2 = ivb2[i];
        ASSERT_EQ((v1!=v2), (ivb1[i]!=ivb2[i]));
        ASSERT_EQ((v1==v2), (ivb1[i]==ivb2[i]));
        ASSERT_EQ((v1<=v2), (ivb1[i]<=ivb2[i]));
        ASSERT_EQ((v1<v2), (ivb1[i]<ivb2[i]));
        ASSERT_EQ((v1>=v2), (ivb1[i]>=ivb2[i]));
        ASSERT_EQ((v1>v2), (ivb1[i]>ivb2[i]));
    }
}

//! Test assign from int_vector_buffer to int_vector_buffer and compare operators
TEST_F(IntVectorBufferTest, Compare)
{
    for (size_type width=1; width<=64; ++width) { // for each possible width
        compare< sdsl::int_vector_buffer<> >(width);
    }
    compare< sdsl::int_vector_buffer<1> >();
    compare< sdsl::int_vector_buffer<8> >();
    compare< sdsl::int_vector_buffer<13> >();
    compare< sdsl::int_vector_buffer<16> >();
    compare< sdsl::int_vector_buffer<32> >();
    compare< sdsl::int_vector_buffer<64> >();
}

template<class t_T>
void test_push_back(size_type width=1)
{
    std::mt19937_64 rng(13);
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;
    {
        t_T ivb(file_name, false, buffersize, width);
        for (size_type i=0; i < 100000; ++i) {
            value_type x = rng() & sdsl::bits::lo_set[ivb.width()];
            ivb.push_back(x);
            ASSERT_EQ(x, (size_type)ivb[i]);
        }
        ivb.close();
    }
    rng.seed(13); // To get the same values
    {
        t_T ivb(file_name, true, buffersize, width);
        for (size_type i=0; i < ivb.size(); ++i) {
            value_type x = rng() & sdsl::bits::lo_set[ivb.width()];
            ASSERT_EQ(x, (size_type)ivb[i]);
        }
    }
}

//! Test push_back
TEST_F(IntVectorBufferTest, PushBack)
{
    for (size_type width=1; width <= 64; ++width) { // for each possible width
        test_push_back< sdsl::int_vector_buffer<> >(width);
    }
    test_push_back< sdsl::int_vector_buffer<1> >();
    test_push_back< sdsl::int_vector_buffer<8> >();
    test_push_back< sdsl::int_vector_buffer<13> >();
    test_push_back< sdsl::int_vector_buffer<16> >();
    test_push_back< sdsl::int_vector_buffer<32> >();
    test_push_back< sdsl::int_vector_buffer<64> >();
}

template<class t_T>
void test_swap(size_type exp_w_1, size_type exp_w_2, std::vector<size_type>& vec_sizes)
{
    std::string file_name = "tmp/int_vector_buffer";
    std::string file_name_tmp = "tmp/int_vector_buffer_tmp";
    for (auto size : vec_sizes) {
        if (size < 1000) {
            t_T ivb(file_name, false, 100, exp_w_1, true);
            for (size_type j=0; j < size; ++j) {
                ivb[j] = j & sdsl::bits::lo_set[exp_w_1];
            }
            ASSERT_EQ(file_name, ivb.filename());
            ASSERT_EQ(size, ivb.size());
            ASSERT_EQ(exp_w_1, ivb.width());
            ASSERT_LT((ivb.buffersize()-100), ivb.width());
            ASSERT_EQ(true, ivb.persistence());
            for (size_type j=0; j < size; ++j) {
                ASSERT_EQ(j & sdsl::bits::lo_set[exp_w_1], (size_type)ivb[j]);
            }
            size_type buffersize_ivb = ivb.buffersize();
            t_T tmp(file_name_tmp, false, 80, 2, false);
            ASSERT_EQ(file_name_tmp, tmp.filename());
            ASSERT_EQ((size_type)0, tmp.size());
            ASSERT_EQ(exp_w_2, tmp.width());
            ASSERT_LT((tmp.buffersize()-80), tmp.width());
            ASSERT_EQ(false, tmp.persistence());
            size_type buffersize_tmp = tmp.buffersize();
            tmp.swap(ivb);
            ASSERT_EQ(file_name_tmp, ivb.filename());
            ASSERT_EQ(file_name, tmp.filename());
            ASSERT_EQ((size_type)0, ivb.size());
            ASSERT_EQ(size, tmp.size());
            ASSERT_EQ(buffersize_tmp, ivb.buffersize());
            ASSERT_EQ(buffersize_ivb, tmp.buffersize());
            ASSERT_EQ(exp_w_2, ivb.width());
            ASSERT_EQ(exp_w_1, tmp.width());
            ASSERT_EQ(false, ivb.persistence());
            ASSERT_EQ(true, tmp.persistence());
            for (size_type j=0; j < tmp.size(); ++j) {
                ASSERT_EQ(j & sdsl::bits::lo_set[exp_w_1], (size_type)tmp[j]);
            }
            tmp.persistence(false);
        }
    }
}

//! Test swap
TEST_F(IntVectorBufferTest, Swap)
{
    for (size_type width=1; width <= 64; ++width) { // for each possible width
        test_swap< sdsl::int_vector_buffer<> >(width, 2, vec_sizes);
    }
    test_swap< sdsl::int_vector_buffer<1> >(1, 1, vec_sizes);
    test_swap< sdsl::int_vector_buffer<8> >(8, 8, vec_sizes);
    test_swap< sdsl::int_vector_buffer<13> >(13, 13, vec_sizes);
    test_swap< sdsl::int_vector_buffer<16> >(16, 16, vec_sizes);
    test_swap< sdsl::int_vector_buffer<32> >(32, 32, vec_sizes);
    test_swap< sdsl::int_vector_buffer<64> >(64, 64, vec_sizes);
}

template<class t_T, class t_V>
void test_file_handling(size_type exp_w)
{
    std::mt19937_64 rng(13);
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;

    // write int_vector-file using int_vector_buffer
    t_T ivb(file_name, false, buffersize, exp_w, true);
    for (size_type i=0; i<100000; ++i) {
        ivb[i] = rng() & sdsl::bits::lo_set[exp_w];
    }
    ASSERT_EQ((size_type)100000, ivb.size());
    ASSERT_EQ(exp_w, ivb.width());
    ivb.close();

    // load int_vector from file and compare it with int_vector_buffer
    t_V iv;
    sdsl::load_from_file(iv, file_name);
    ASSERT_EQ((size_type)100000, iv.size());
    t_T ivb2(file_name, true, buffersize, 5, false);
    ASSERT_EQ((size_type)100000, ivb2.size());
    ASSERT_EQ(exp_w, ivb2.width());
    rng.seed(13);
    for (size_type i=0; i<iv.size(); ++i) {
        size_type exp_val = rng() & sdsl::bits::lo_set[exp_w];
        ASSERT_EQ(exp_val, (size_type)ivb2[i]);
        ASSERT_EQ(exp_val, (size_type)iv[i]);
    }
    ivb2.close();
}

template<class t_T>
void test_plain_file_handling(uint8_t exp_w)
{
    std::mt19937_64 rng(13);
    std::string file_name = "tmp/int_array";
    size_type buffersize = 1024;

    // write plain array using int_vector_buffer
    t_T ivb(file_name, false, buffersize, exp_w, true, true);
    for (size_type i=0; i<100000; ++i) {
        ivb[i] = rng();
    }
    ivb.close();

    // load plain_array and compare it with int_vector_buffer
    sdsl::int_vector<> iv;
    sdsl::load_vector_from_file(iv, file_name, exp_w/8);
    ASSERT_EQ((size_type)100000, iv.size());        // Check size
    ASSERT_EQ(exp_w, iv.width());        // Check width
    t_T ivb2(file_name, true, buffersize, exp_w, false, true); // Load plain data with int_vector_buffer
    ASSERT_EQ((size_type)100000, ivb2.size());      // Check size
    ASSERT_EQ(exp_w, ivb2.width());      // Check width
    rng.seed(13);
    for (size_type i=0; i<iv.size(); ++i) {                    // Check content
        size_type exp_val = rng() & sdsl::bits::lo_set[exp_w];
        ASSERT_EQ(exp_val, (size_type)ivb2[i]);
        ASSERT_EQ(exp_val, iv[i]);
    }
    ivb2.close();
}

//! Test the persistence in constructor, opening of existing file, int_vector-format and close()
TEST_F(IntVectorBufferTest, FileHandling)
{
    for (size_type width=1; width <= 64; ++width) { // for each possible width
        test_file_handling< sdsl::int_vector_buffer<>, sdsl::int_vector<> >(width);
    }
    test_file_handling< sdsl::int_vector_buffer<>, sdsl::int_vector<> >(42);
    test_file_handling< sdsl::int_vector_buffer<1>, sdsl::int_vector<1> >(1);
    test_file_handling< sdsl::int_vector_buffer<8>, sdsl::int_vector<8> >(8);
    test_file_handling< sdsl::int_vector_buffer<13>, sdsl::int_vector<13> >(13);
    test_file_handling< sdsl::int_vector_buffer<16>, sdsl::int_vector<16> >(16);
    test_file_handling< sdsl::int_vector_buffer<32>, sdsl::int_vector<32> >(32);
    test_file_handling< sdsl::int_vector_buffer<64>, sdsl::int_vector<64> >(64);
    test_plain_file_handling< sdsl::int_vector_buffer<> >(8);
    test_plain_file_handling< sdsl::int_vector_buffer<> >(16);
    test_plain_file_handling< sdsl::int_vector_buffer<> >(32);
    test_plain_file_handling< sdsl::int_vector_buffer<> >(64);
    test_plain_file_handling< sdsl::int_vector_buffer<8> >(8);
    test_plain_file_handling< sdsl::int_vector_buffer<16> >(16);
    test_plain_file_handling< sdsl::int_vector_buffer<32> >(32);
    test_plain_file_handling< sdsl::int_vector_buffer<64> >(64);
}

template<class t_T>
void test_random_access(size_type width=1)
{
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, 999);
    auto dice = bind(distribution, rng);
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 100;
    t_T ivb(file_name, false, buffersize, width);
    ASSERT_EQ((size_type)0, ivb.size());
    ivb[999] = 999 & sdsl::bits::lo_set[ivb.width()];
    ASSERT_EQ(999 & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[999]);
    ASSERT_EQ((size_type)1000, ivb.size());
    for (size_type i=0; i < 999; ++i) {
        ASSERT_EQ((size_type)0, (size_type)ivb[i]);
    }
    // Test random access write
    for (size_type i=0; i < 1000; ++i) {
        value_type x = dice();
        ivb[x] = x & sdsl::bits::lo_set[ivb.width()];
    }
    ASSERT_EQ((size_type)1000, ivb.size());
    for (size_type i=0; i < ivb.size(); ++i) {
        if (ivb[i]) {
            ASSERT_EQ(i & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        }
    }
    // Test random access with different buffersize
    buffersize = 50;
    ivb.buffersize(buffersize);
    ASSERT_EQ((size_type)1000, ivb.size());
    for (size_type i=0; i < 1000; ++i) {
        value_type x = dice();
        ivb[x] = x & sdsl::bits::lo_set[ivb.width()];
    }
    for (size_type i=0; i < ivb.size(); ++i) {
        if (ivb[i]) {
            ASSERT_EQ(i & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        }
    }
    // Test random access read
    for (size_type i=0; i < ivb.size(); ++i) {
        value_type idx = dice();
        if (ivb[idx]) {
            ASSERT_EQ(idx & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[idx]);
        }
    }
}

//! Test RandomAcces, which should not be done in practice because it is expected to be very slow
TEST_F(IntVectorBufferTest, RandomAccess)
{
    for (size_type width=1; width <= 64; ++width) {
        test_random_access< sdsl::int_vector_buffer<> >(width);
    }
    test_random_access< sdsl::int_vector_buffer<1> >();
    test_random_access< sdsl::int_vector_buffer<8> >();
    test_random_access< sdsl::int_vector_buffer<16> >();
    test_random_access< sdsl::int_vector_buffer<32> >();
    test_random_access< sdsl::int_vector_buffer<64> >();
}

template<class t_T>
void test_move(size_type constructor_width)
{
    std::string file_name = "tmp/int_vector_buffer";
    uint64_t value = 3;

    {
        // Test MoveConstructor
        std::vector<t_T> v;
        for (uint64_t i=0; i<4; i++) {
            v.push_back(t_T(file_name+sdsl::util::to_string(i), false, 1000, i+1, false));
            v[i][5] = value;
            ASSERT_EQ(value & sdsl::bits::lo_set[v[i].width()], (size_type)v[i][5]);
        }
        for (uint64_t i=0; i<v.size(); i++) {
            v[i][5] = value;
            ASSERT_EQ(value & sdsl::bits::lo_set[v[i].width()], (size_type)v[i][5]);
        }
        {
            t_T tmp(file_name, false, 1000, 13, false);
            ASSERT_EQ((size_type)0, tmp.size());
            ASSERT_LT((tmp.buffersize()-(1000)), tmp.width());
            size_type buffersize = tmp.buffersize();
            t_T ivb((t_T&&)tmp);
            ASSERT_EQ(file_name, ivb.filename());
            ASSERT_EQ((size_type)0, ivb.size());
            ASSERT_EQ(buffersize, ivb.buffersize());
            ASSERT_EQ("", tmp.filename());
            ASSERT_EQ((size_type)0, tmp.size());
            ASSERT_EQ((size_type)0, tmp.buffersize());
            ASSERT_EQ(constructor_width, tmp.width());
        }
    }

    {
        // Test MoveAssignment
        std::vector<t_T> v(4);
        for (uint64_t i=0; i<v.size(); i++) {
            v[i] = t_T(file_name+sdsl::util::to_string(i), false, 1000, i+1, false);
            v[i][5] = value;
            ASSERT_EQ(value & sdsl::bits::lo_set[v[i].width()], (size_type)v[i][5]);
        }
        for (uint64_t i=0; i<v.size(); i++) {
            t_T tmp(file_name+sdsl::util::to_string(i), false, 1000, i+1, false);
            v[i] = (t_T&&)tmp;
            ASSERT_EQ("", tmp.filename());
            ASSERT_EQ((size_type)0, tmp.size());
            ASSERT_EQ((size_type)0, tmp.buffersize());
            ASSERT_EQ(constructor_width, tmp.width());
            bool persistence = false;
            ASSERT_EQ(persistence, (bool)tmp.persistence());
            ASSERT_EQ(file_name+sdsl::util::to_string(i), v[i].filename());
        }
    }

}

//! Test MoveConstructor and MoveAssignment
TEST_F(IntVectorBufferTest, Move)
{
    test_move< sdsl::int_vector_buffer<> >(64);
    test_move< sdsl::int_vector_buffer<1> >(1);
    test_move< sdsl::int_vector_buffer<8> >(8);
    test_move< sdsl::int_vector_buffer<16> >(16);
    test_move< sdsl::int_vector_buffer<32> >(32);
    test_move< sdsl::int_vector_buffer<64> >(64);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
