#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::int_vector<>::value_type value_type;

std::string temp_dir;

// The fixture for testing class int_vector.
class int_vector_buffer_test : public ::testing::Test
{
    protected:

        int_vector_buffer_test() {}

        virtual ~int_vector_buffer_test() {}

        virtual void SetUp()
        {
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
    static_assert(std::is_default_constructible<t_T>::value, "Type is not default constructible");
    static_assert(std::is_move_constructible<t_T>::value, "Type is not move constructible");
    static_assert(std::is_move_assignable<t_T>::value, "Type is not move assignable");
    std::string file_name = temp_dir+"/int_vector_buffer";
    {
        // Default constructor
        t_T ivb;
        ASSERT_FALSE(ivb.is_open());                      // int_vector_buffer is not ready for IO since no filename was given
        ASSERT_EQ("", ivb.filename());                    // filename should be empty
        ASSERT_EQ((size_type)0, ivb.size());              // size should be 0
        ASSERT_EQ((uint8_t)template_width, ivb.width());  // default width of each element should be template_width bits
    }
    {
        // Constructor with 2 Parameters
        t_T ivb(file_name, std::ios::out);
        ASSERT_TRUE(ivb.is_open());                             // int_vector_buffer should be open
        ASSERT_EQ(file_name, ivb.filename());                   // filename should be file_name
        ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
        ASSERT_EQ((uint8_t)template_width, ivb.width());        // default width of each element should be template_width bits
        ivb.close(true);
    }
    {
        // Constructor with 4 Parameters
        size_type buffersize = 1024*1024;
        t_T ivb(file_name, std::ios::out, buffersize, constructor_width);
        ASSERT_TRUE(ivb.is_open());                             // int_vector_buffer should be open
        ASSERT_EQ(file_name, ivb.filename());                   // filename should be file_name
        ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
        ASSERT_EQ((uint8_t)exp_width, ivb.width());             // default width of each element should be i bits
        ivb.close(true);
    }
}

//! Test constructors
TEST_F(int_vector_buffer_test, constructors)
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
void test_assign_and_modify(size_type width=1)
{
    std::mt19937_64 rng(13), rng2;
    std::string file_name = temp_dir+"/int_vector_buffer";
    size_type buffersize = 1024;
    t_T ivb(file_name, std::ios::out, buffersize, width);
    for (size_type i=0; i < 100000; ++i) {
        value_type exp_v = rng();
        ivb[i] = exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Assign Operator failed";
    }
    rng.seed(13); // To get the same values
    for (size_type i=0; i < 100000; ++i) {
        value_type exp_v = rng(), tmp = rng2();
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Assign Operator failed";
        ivb[i] += tmp;
        exp_v += tmp;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Add Assign Operator failed";
        ivb[i] -= tmp;
        exp_v -= tmp;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Subtract Assign Operator failed";
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]++) << "Postfix Increment Operator failed";
        exp_v++;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Postfix Increment Operator failed";
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]--) << "Postfix Decrement Operator failed";
        exp_v--;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Postfix Decrement Operator failed";
        ++exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)++ivb[i]) << "Prefix Increment Operator failed";
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Prefix Increment Operator failed";
        --exp_v;
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)--ivb[i]) << "Prefix Decrement Operator failed";
        ASSERT_EQ(exp_v & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]) << "Prefix Decrement Operator failed";
    }
    ivb.close(true);
}

template<>
void test_assign_and_modify<sdsl::int_vector_buffer<1>>(size_type width)
{
    std::mt19937_64 rng(13);
    std::uniform_int_distribution<uint64_t> distribution(0, 9);

    std::string file_name = temp_dir+"/int_vector_buffer";
    size_type buffersize = 1024;
    sdsl::int_vector_buffer<1> ivb(file_name, std::ios::out, buffersize, width);
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
    ivb.close(true);
}

TEST_F(int_vector_buffer_test, assign_and_modify_element)
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
    std::string file_name_1 = temp_dir+"/int_vector_buffer_1";
    std::string file_name_2 = temp_dir+"/int_vector_buffer_2";
    size_type buffersize = 1024;
    t_T ivb1(file_name_1, std::ios::out, buffersize, width);
    t_T ivb2(file_name_2, std::ios::out, buffersize, width);
    for (size_type i=0; i < 1000; ++i) {
        ivb1[i] = rng() & sdsl::bits::lo_set[ivb1.width()];
        ivb2[i] = ivb1[i];
        ASSERT_TRUE(ivb1[i]==ivb2[i]);
        ASSERT_FALSE(ivb1[i]!=ivb2[i]);
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
    ivb1.close(true);
    ivb2.close(true);
}

//! Test assign from int_vector_buffer to int_vector_buffer and compare operators
TEST_F(int_vector_buffer_test, compare)
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
void test_sequential_access(size_type width=1)
{
    std::mt19937_64 rng;
    std::string file_name = temp_dir+"/int_vector_buffer";
    size_type buffersize = 1024;
    size_type size = 100000;
    // fill ivb with push_back()
    {
        rng.seed(13); // To get the same values
        t_T ivb(file_name, std::ios::out, buffersize, width);
        for (size_type i=0; i < size; ++i) {
            value_type x = rng() & sdsl::bits::lo_set[ivb.width()];
            ivb.push_back(x);
            ASSERT_EQ(x, (size_type)ivb[i]);
        }
    }
    // verify values
    {
        rng.seed(13); // To get the same values
        t_T ivb(file_name, std::ios::in, buffersize, width);
        ASSERT_EQ(size, ivb.size());
        for (size_type i=0; i < ivb.size(); ++i) {
            value_type x = rng() & sdsl::bits::lo_set[ivb.width()];
            ASSERT_EQ(x, (size_type)ivb[i]);
        }
    }
    // iterate over ivb, verify and change values
    {
        rng.seed(13); // To get the same values
        t_T ivb(file_name, std::ios::in, buffersize, width);
        for (typename t_T::iterator it = ivb.begin(); it != ivb.end(); ++it) {
            value_type x = rng() & sdsl::bits::lo_set[ivb.width()];
            ASSERT_EQ(x, (size_type)*it);
            *it = (x+1) & sdsl::bits::lo_set[ivb.width()];
        }
    }
    // iterator test
    {
        t_T ivb(file_name, std::ios::in, buffersize, width);
        typename t_T::iterator it = ivb.begin();
        ASSERT_EQ(ivb[0], *it++);
        ASSERT_EQ(ivb[1], *it);
        it += 1;
        ASSERT_EQ(ivb[2], *it);
        it -= 1;
        ASSERT_EQ(ivb[1], *it);
        it -= -1;
        ASSERT_EQ(ivb[2], *it);
        it += -1;
        ASSERT_EQ(ivb[1], *it);
        ASSERT_EQ(ivb[2], *(++it));
        typename t_T::iterator it2 = ivb.end();
        --it2;
        ASSERT_EQ(ivb[ivb.size()-1], *it2--);
        ASSERT_EQ(ivb[ivb.size()-2], *it2);
        ASSERT_EQ(ivb[ivb.size()-3], *(--it2));
    }

    // verify changed values
    {
        rng.seed(13); // To get the same values
        t_T ivb(file_name, std::ios::in, buffersize, width);
        ASSERT_EQ(size, ivb.size());
        for (size_type i=0; i < ivb.size(); ++i) {
            value_type x = (rng()+1) & sdsl::bits::lo_set[ivb.width()];
            ASSERT_EQ(x, (size_type)ivb[i]) << "???";
        }
        ivb.close(true);
    }
}

//! Test SequentialAccess: push_back and iterators
TEST_F(int_vector_buffer_test, sequential_access)
{
    for (size_type width=1; width <= 64; ++width) { // for each possible width
        test_sequential_access< sdsl::int_vector_buffer<> >(width);
    }
    test_sequential_access< sdsl::int_vector_buffer<1> >();
    test_sequential_access< sdsl::int_vector_buffer<8> >();
    test_sequential_access< sdsl::int_vector_buffer<13> >();
    test_sequential_access< sdsl::int_vector_buffer<16> >();
    test_sequential_access< sdsl::int_vector_buffer<32> >();
    test_sequential_access< sdsl::int_vector_buffer<64> >();
}


template<class t_T>
void test_random_access(size_type width=1)
{
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, 999);
    auto dice = bind(distribution, rng);
    std::string file_name = temp_dir+"/int_vector_buffer";
    size_type buffersize = 100;
    t_T ivb(file_name, std::ios::out, buffersize, width);
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
    ivb.close(true);
}

//! Test RandomAcces, which should not be done in practice because it is expected to be very slow
TEST_F(int_vector_buffer_test, random_access)
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


template<class t_T, class t_V>
void test_file_handling(size_type exp_w)
{
    std::mt19937_64 rng(13);
    std::string file_name = temp_dir+"/int_vector_buffer";
    size_type buffersize = 1024;

    // write an int_vector-file using int_vector_buffer
    t_T ivb(file_name, std::ios::out, buffersize, exp_w);
    for (size_type i=0; i<100000; ++i) {
        ivb[i] = rng() & sdsl::bits::lo_set[exp_w];
    }
    ivb.close();

    // load int_vector from file and compare it with int_vector_buffer
    t_V iv;
    sdsl::load_from_file(iv, file_name);
    ASSERT_EQ((size_type)100000, iv.size());
    t_T ivb2(file_name, std::ios::in, buffersize, 5);
    ASSERT_EQ((size_type)100000, ivb2.size());
    ASSERT_EQ(exp_w, ivb2.width());
    rng.seed(13);
    for (size_type i=0; i<iv.size(); ++i) {
        size_type exp_val = rng() & sdsl::bits::lo_set[exp_w];
        ASSERT_EQ(exp_val, (size_type)ivb2[i]);
        ASSERT_EQ(exp_val, (size_type)iv[i]);
    }
    ivb2.close(true);
}

template<class t_T>
void test_plain_file_handling(uint8_t exp_w)
{
    std::mt19937_64 rng(13);
    std::string file_name = temp_dir+"/int_array";
    size_type buffersize = 1024;

    // write plain array using int_vector_buffer
    t_T ivb(file_name, std::ios::out, buffersize, exp_w, true);
    for (size_type i=0; i<100000; ++i) {
        ivb[i] = rng();
    }
    ivb.close();

    // load plain_array and compare it with int_vector_buffer
    sdsl::int_vector<> iv;
    sdsl::load_vector_from_file(iv, file_name, exp_w/8);
    ASSERT_EQ((size_type)100000, iv.size()) << "written plain int-array has wrong size";
    ASSERT_EQ(exp_w, iv.width()) << "written plain int-array has wrong width";
    t_T ivb2(file_name, std::ios::in, buffersize, exp_w, true); // Load plain data with int_vector_buffer
    ASSERT_EQ((size_type)100000, ivb2.size());                  // Check size
    ASSERT_EQ(exp_w, ivb2.width());                             // Check width
    rng.seed(13);
    for (size_type i=0; i<iv.size(); ++i) {                     // Check content
        size_type exp_val = rng() & sdsl::bits::lo_set[exp_w];
        ASSERT_EQ(exp_val, (size_type)ivb2[i]);
        ASSERT_EQ(exp_val, iv[i]);
    }
    ivb2.close(true);
}

//! Test opening of existing file, int_vector-format and close()
TEST_F(int_vector_buffer_test, file_handling)
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
void test_swap(size_type exp_w_ivb1, size_type exp_w_ivb2, std::vector<size_type>& vec_sizes)
{
    std::string file_name_1 = temp_dir+"/int_vector_buffer_1";
    std::string file_name_2 = temp_dir+"/int_vector_buffer_2";
    for (auto size : vec_sizes) {
        if (size < 1000) {
            // Create, fill and verify ivb1
            t_T ivb1(file_name_1, std::ios::out, 100, exp_w_ivb1);
            for (size_type j=0; j < size; ++j) {
                ivb1[j] = j & sdsl::bits::lo_set[exp_w_ivb1];
            }
            ASSERT_TRUE(ivb1.is_open());
            ASSERT_EQ(file_name_1, ivb1.filename());
            ASSERT_EQ(size, ivb1.size());
            ASSERT_EQ(exp_w_ivb1, ivb1.width());
            for (size_type j=0; j < size; ++j) {
                ASSERT_EQ(j & sdsl::bits::lo_set[exp_w_ivb1], (size_type)ivb1[j]);
            }
            size_type buffersize_ivb1 = ivb1.buffersize();

            // Create and verify ivb2
            t_T ivb2(file_name_2, std::ios::out, 80, 2);
            ASSERT_TRUE(ivb2.is_open());
            ASSERT_EQ(file_name_2, ivb2.filename());
            ASSERT_EQ((size_type)0, ivb2.size());
            ASSERT_EQ(exp_w_ivb2, ivb2.width());
            size_type buffersize_ivb2 = ivb2.buffersize();

            // Swap ivb1 and ivb2
            ivb1.swap(ivb2);

            // Check ivb1
            ASSERT_TRUE(ivb1.is_open());
            ASSERT_EQ(file_name_2, ivb1.filename());
            ASSERT_EQ((size_type)0, ivb1.size());
            ASSERT_EQ(buffersize_ivb2, ivb1.buffersize());
            ASSERT_EQ(exp_w_ivb2, ivb1.width());
            ivb1.close(true);

            // Check ivb2
            ASSERT_TRUE(ivb2.is_open());
            ASSERT_EQ(file_name_1, ivb2.filename());
            ASSERT_EQ(size, ivb2.size());
            ASSERT_EQ(buffersize_ivb1, ivb2.buffersize());
            ASSERT_EQ(exp_w_ivb1, ivb2.width());
            for (size_type j=0; j < ivb2.size(); ++j) {
                ASSERT_EQ(j & sdsl::bits::lo_set[exp_w_ivb1], (size_type)ivb2[j]);
            }
            ivb2.close(true);
        }
    }
}

//! Test swap
TEST_F(int_vector_buffer_test, swap)
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


template<class t_T>
void test_move(size_type constructor_width)
{
    std::string file_name = temp_dir+"/int_vector_buffer";
    std::mt19937_64 rng;
    size_type numbers = 10000;
    // Test MoveConstructor in vector
    {
        std::vector<t_T> v;
        for (uint64_t i=0; i<4; i++) {
            v.push_back(t_T(file_name+sdsl::util::to_string(i), std::ios::out, 1000, i+1));
            ASSERT_TRUE(v[i].is_open());
            // fill
            rng.seed(13);
            for (uint64_t j=0; j<numbers; ++j) {
                size_type value = rng() & sdsl::bits::lo_set[v[i].width()];
                v[i][j] = value;
                ASSERT_EQ(value, (size_type)v[i][j]);
            }
        }
        // Check if all values are correct
        for (uint64_t i=0; i<v.size(); i++) {
            rng.seed(13);
            for (uint64_t j=0; j<numbers; ++j) {
                size_type value = rng() & sdsl::bits::lo_set[v[i].width()];
                ASSERT_EQ(value, (size_type)v[i][j]);
            }
            v[i].close(true);
        }
    }
    // Test MoveConstructor
    {
        t_T tmp(file_name, std::ios::out, 1000, 13);
        rng.seed(13);
        for (uint64_t j=0; j<numbers; ++j) {
            size_type value = rng() & sdsl::bits::lo_set[tmp.width()];
            tmp[j] = value;
        }
        // MoveConstructor
        t_T ivb(std::move(tmp));

        // Check ivb
        ASSERT_TRUE(ivb.is_open());
        rng.seed(13);
        for (uint64_t j=0; j<numbers; ++j) {
            size_type value = rng() & sdsl::bits::lo_set[ivb.width()];
            ASSERT_EQ(value, (size_type)ivb[j]);
        }
        ASSERT_EQ(file_name, ivb.filename());
        ASSERT_EQ((size_type)numbers, ivb.size());
        ivb.close(true);

        // Check tmp
        ASSERT_FALSE(tmp.is_open());
        ASSERT_EQ("", tmp.filename());
        ASSERT_EQ((size_type)0, tmp.size());
        ASSERT_EQ(constructor_width, tmp.width());
    }

    // Test MoveAssignment
    {
        std::vector<t_T> v(4);
        for (uint64_t i=0; i<v.size(); i++) {
            v[i] = t_T(file_name+sdsl::util::to_string(i), std::ios::out, 1000, i+1);
            ASSERT_TRUE(v[i].is_open());
            rng.seed(13);
            for (uint64_t j=0; j<numbers; ++j) {
                size_type value = rng() & sdsl::bits::lo_set[v[i].width()];
                v[i][j] = value;
            }
            rng.seed(13);
            for (uint64_t j=0; j<numbers; ++j) {
                size_type value = rng() & sdsl::bits::lo_set[v[i].width()];
                ASSERT_EQ(value, (size_type)v[i][j]);
            }
            v[i].close(true);
        }
        for (uint64_t i=0; i<v.size(); i++) {
            t_T tmp(file_name+sdsl::util::to_string(i), std::ios::out, 1000, i+1);
            size_type buffersize = tmp.buffersize();
            rng.seed(13);
            for (uint64_t j=0; j<numbers; ++j) {
                size_type value = rng() & sdsl::bits::lo_set[v[i].width()];
                tmp[j] = value;
            }
            v[i] = std::move(tmp);

            // Check v[i]
            ASSERT_TRUE(v[i].is_open());
            rng.seed(13);
            for (uint64_t j=0; j<numbers; ++j) {
                size_type value = rng() & sdsl::bits::lo_set[v[i].width()];
                ASSERT_EQ(value, (size_type)v[i][j]);
            }
            ASSERT_EQ(file_name+sdsl::util::to_string(i), v[i].filename());
            ASSERT_EQ((size_type)numbers, v[i].size());
            ASSERT_EQ(buffersize, v[i].buffersize());
            v[i].close(true);

            // Check tmp
            ASSERT_FALSE(tmp.is_open());
            ASSERT_EQ("", tmp.filename());
            ASSERT_EQ((size_type)0, tmp.size());
            ASSERT_EQ(constructor_width, tmp.width());
        }
    }
}

//! Test MoveConstructor and MoveAssignment
TEST_F(int_vector_buffer_test, move)
{
    test_move< sdsl::int_vector_buffer<> >(64);
    test_move< sdsl::int_vector_buffer<1> >(1);
    test_move< sdsl::int_vector_buffer<8> >(8);
    test_move< sdsl::int_vector_buffer<16> >(16);
    test_move< sdsl::int_vector_buffer<32> >(32);
    test_move< sdsl::int_vector_buffer<64> >(64);
}


template<class t_T>
void test_reset(std::vector<size_type>& vec_sizes, size_type width=1)
{
    std::mt19937_64 rng(13);
    std::string file_name = temp_dir+"/int_vector_buffer";
    size_type buffersize = 1024;
    for (auto size : vec_sizes) {
        if (size < 1000) {
            t_T ivb(file_name, std::ios::out, buffersize, width);
            for (size_type j=0; j < size; ++j) {
                ivb[j] = rng();
            }
            ASSERT_EQ(file_name, ivb.filename());
            ASSERT_EQ(size, ivb.size());
            size_type bsize = ivb.buffersize();
            ivb.reset();                           // reset should delete all content
            ASSERT_EQ(file_name, ivb.filename());  // same filename as before
            ASSERT_EQ((size_type)0, ivb.size());   // all content removed
            ASSERT_EQ(bsize, ivb.buffersize());    // same buffersize as before
            {
                std::ifstream ifile(file_name, std::ios::in|std::ios::binary|std::ios::ate);
                size_type file_end = ifile.tellg();
                ifile.close();
                ASSERT_EQ((size_type)0, file_end);  // size of file after reset is 0
            }
            ivb.close(true);
        }
    }
}

//! Test reset
TEST_F(int_vector_buffer_test, Reset)
{
    for (size_type width=1; width <= 64; ++width) {
        test_reset< sdsl::int_vector_buffer<> >(vec_sizes, width);
    }
    test_reset< sdsl::int_vector_buffer<1> >(vec_sizes);
    test_reset< sdsl::int_vector_buffer<8> >(vec_sizes);
    test_reset< sdsl::int_vector_buffer<16> >(vec_sizes);
    test_reset< sdsl::int_vector_buffer<32> >(vec_sizes);
    test_reset< sdsl::int_vector_buffer<64> >(vec_sizes);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        std::cout << "Usage: " << argv[0] << " tmp_dir" << std::endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    temp_dir = argv[1];
    return RUN_ALL_TESTS();
}
