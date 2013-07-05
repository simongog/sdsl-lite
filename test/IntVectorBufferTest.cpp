#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector_buffer.hpp"
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
class IntVectorBufferTest : public ::testing::Test
{
    protected:

        IntVectorBufferTest() {}

        virtual ~IntVectorBufferTest() {}

        virtual void SetUp() {
            for (size_type i=0; i<128; ++i) {
                vec_sizes.push_back(rand() % 100000);
            }
            for (size_type i=0; i < 10; ++i) {
                vec_sizes.push_back(rand()%10000000);
            }
        }

        virtual void TearDown() {}

        std::vector<size_type> vec_sizes = {0, 64, 65, 127, 128}; // different sizes for the vectors
};

template<class t_T>
void test_default_constructor(size_type exp_w)
{
    {
        t_T ivb;
        ASSERT_EQ((size_type)0, ivb.size());         // size should be 0
        ASSERT_EQ((size_type)0, ivb.buffersize());   // buffersize should be 0
        ASSERT_EQ((uint8_t)exp_w, ivb.width());      // default width of each element should be exp_w bits
    }
}

//! Test the default constructor
TEST_F(IntVectorBufferTest, DefaultConstruct)
{
    test_default_constructor< sdsl::int_vector_buffer<> >(64);
    test_default_constructor< sdsl::int_vector_buffer<13> >(13);
    test_default_constructor< sdsl::int_vector_buffer<1> >(1);
    test_default_constructor< sdsl::int_vector_buffer<8> >(8);
    test_default_constructor< sdsl::int_vector_buffer<16> >(16);
    test_default_constructor< sdsl::int_vector_buffer<32> >(32);
    test_default_constructor< sdsl::int_vector_buffer<64> >(64);
}

template<class t_T>
void test_two_param_constructor(size_type exp_w)
{
    std::string file_name = "tmp/int_vector_buffer";
    {
        t_T ivb(file_name, false);
        ASSERT_EQ((size_type)0, ivb.size());                   // size should be 0
        ASSERT_EQ((uint8_t)exp_w, ivb.width());                // default width of each element should be exp_w bits
        ASSERT_LT((ivb.buffersize()-(size_type)(1024*1024)), ivb.width());  // actual buffersize() is not more than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
    }
}

//! Test the parametrized constructor
TEST_F(IntVectorBufferTest, ConstructorTwoParam)
{
    test_two_param_constructor< sdsl::int_vector_buffer<> >(64);
    test_two_param_constructor< sdsl::int_vector_buffer<13> >(13);
    test_two_param_constructor< sdsl::int_vector_buffer<1> >(1);
    test_two_param_constructor< sdsl::int_vector_buffer<8> >(8);
    test_two_param_constructor< sdsl::int_vector_buffer<16> >(16);
    test_two_param_constructor< sdsl::int_vector_buffer<32> >(32);
    test_two_param_constructor< sdsl::int_vector_buffer<64> >(64);
}

template<class t_T>
void test_buffersize(std::vector<size_type>& vec_sizes)
{
    std::string file_name = "tmp/int_vector_buffer";
    // Test constructor with different buffersizes
    for (size_type i=0; i < vec_sizes.size(); ++i) {
        t_T ivb(file_name, false, vec_sizes[i]);
        ASSERT_EQ((size_type)0, ivb.size());
        ASSERT_LT((ivb.buffersize()-vec_sizes[i]), ivb.width());  // actual buffersize() is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
    }
    {
        // Test buffersize() method
        t_T ivb(file_name, false);
        for (size_type i=0; i < vec_sizes.size(); ++i) {
            ivb.buffersize(vec_sizes[i]);
            ASSERT_EQ((size_type)0, ivb.size());
            ASSERT_LT((ivb.buffersize()-vec_sizes[i]), ivb.width());  // actual buffersize() is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
    }
    {
        // Test buffersize() method
        t_T ivb(file_name, false, 0);
        for (size_type i=0; i < vec_sizes.size(); ++i) {
            ivb.buffersize(vec_sizes[i]);
            ASSERT_EQ((size_type)0, ivb.size());
            ASSERT_LT((ivb.buffersize()-vec_sizes[i]), ivb.width());  // actual buffersize() is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
    }
}

//! Test the buffersize in constructor and buffersize()
TEST_F(IntVectorBufferTest, Buffersize)
{
    test_buffersize< sdsl::int_vector_buffer<> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<13> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<1> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<8> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<16> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<32> >(vec_sizes);
    test_buffersize< sdsl::int_vector_buffer<64> >(vec_sizes);
}

//! Test the parametrized constructor with four parameters
TEST_F(IntVectorBufferTest, ConstructorFourParam)
{
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024*1024;
    for (size_type width=1; width<=64; ++width) {
        {
            sdsl::int_vector_buffer<> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)width, ivb.width());                 // default width of each element should be i bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
        {
            sdsl::int_vector_buffer<13> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)13, ivb.width());                    // default width of each element should be i bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
        {
            sdsl::int_vector_buffer<42> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)42, ivb.width());                    // default width of each element should be i bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
        {
            sdsl::int_vector_buffer<1> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)1, ivb.width());                     // default width of each element should be 1 bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
        {
            sdsl::int_vector_buffer<8> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)8, ivb.width());                     // default width of each element should be 8 bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
        {
            sdsl::int_vector_buffer<16> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)16, ivb.width());                    // default width of each element should be 16 bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
        {
            sdsl::int_vector_buffer<32> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)32, ivb.width());                    // default width of each element should be 32 bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
        {
            sdsl::int_vector_buffer<64> ivb(file_name, false, buffersize, width);
            ASSERT_EQ((size_type)0, ivb.size());                    // size should be 0
            ASSERT_EQ((uint8_t)64, ivb.width());                    // default width of each element should be 64 bits
            ASSERT_LT((ivb.buffersize()-buffersize), ivb.width());  // actual buffersize is less than 8 elements bigger than given buffersize ( buffersize()*8 < buffersize*8+8*width() )
        }
    }
}

template<class t_T>
void test_assign(unsigned char w=1)
{
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;
    t_T ivb(file_name, false, buffersize, w);
    for (size_type i=0; i < 100000; ++i) {
        value_type x = rand() & sdsl::bits::lo_set[ivb.width()];
        ivb[i] = x;
        ASSERT_EQ(x, (size_type)ivb[i]);
    }
}

TEST_F(IntVectorBufferTest, AssignElement)
{
    for (unsigned char w=1; w <= 64; ++w) { // for each possible width
        test_assign< sdsl::int_vector_buffer<> >(w);
    }
    test_assign< sdsl::int_vector_buffer<1> >();
    test_assign< sdsl::int_vector_buffer<8> >();
    test_assign< sdsl::int_vector_buffer<16> >();
    test_assign< sdsl::int_vector_buffer<32> >();
    test_assign< sdsl::int_vector_buffer<64> >();
}

template<class t_T>
void test_push_back(unsigned char w=1)
{
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;

    uint64_t seed = time(NULL);
    {
        srand(seed);
        t_T ivb(file_name, false, buffersize, w);
        for (size_type i=0; i < 100000; ++i) {
            value_type x = rand() & sdsl::bits::lo_set[ivb.width()];
            ivb.push_back(x);
            ASSERT_EQ(x, (size_type)ivb[i]);
        }
        ivb.close();
    }
    {
        srand(seed);
        t_T ivb(file_name, true, buffersize, w);
        for (size_type i=0; i < ivb.size(); ++i) {
            value_type x = rand() & sdsl::bits::lo_set[ivb.width()];
            ASSERT_EQ(x, (size_type)ivb[i]);
        }
    }
}

TEST_F(IntVectorBufferTest, PushBack)
{
    for (unsigned char w=1; w <= 64; ++w) { // for each possible width
        test_push_back< sdsl::int_vector_buffer<> >(w);
    }
    test_push_back< sdsl::int_vector_buffer<1> >();
    test_push_back< sdsl::int_vector_buffer<8> >();
    test_push_back< sdsl::int_vector_buffer<16> >();
    test_push_back< sdsl::int_vector_buffer<32> >();
    test_push_back< sdsl::int_vector_buffer<64> >();
}

template<class t_T>
void test_swap(size_type exp_w_1, size_type exp_w_2, std::vector<size_type>& vec_sizes)
{
    std::string file_name = "tmp/int_vector_buffer";
    std::string file_name_tmp = "tmp/int_vector_buffer_tmp";
    for (size_type i=0; i < vec_sizes.size()-10; ++i) {
        t_T ivb(file_name, false, 1000, 11, true);
        ASSERT_EQ((size_type)exp_w_1, ivb.width());
        ASSERT_LT((ivb.buffersize()-(size_type)1000), ivb.width());
        size_type buffersize_ivb = ivb.buffersize();
        ASSERT_EQ(true, ivb.persistence());
        for (size_type j=0; j < vec_sizes[i]; ++j) {
            ivb[j] = j & sdsl::bits::lo_set[exp_w_1];
            ASSERT_EQ(j&sdsl::bits::lo_set[exp_w_1], (size_type)ivb[j]);
        }
        ASSERT_EQ(vec_sizes[i], ivb.size());
        t_T tmp(file_name_tmp, false, 800, 2, false);
        ASSERT_EQ((size_type)exp_w_2, tmp.width());
        ASSERT_LT((tmp.buffersize()-(size_type)800), tmp.width());
        size_type buffersize_tmp = tmp.buffersize();
        ASSERT_EQ(false, tmp.persistence());
        ASSERT_EQ((size_type)0, tmp.size());
        tmp.swap(ivb);
        ASSERT_EQ((size_type)0, ivb.size());
        ASSERT_EQ(vec_sizes[i], tmp.size());
        ASSERT_EQ(buffersize_tmp, ivb.buffersize());
        ASSERT_EQ(buffersize_ivb, tmp.buffersize());
        ASSERT_EQ((size_type)exp_w_2, ivb.width());
        ASSERT_EQ((size_type)exp_w_1, tmp.width());
        ASSERT_EQ(false, ivb.persistence());
        ASSERT_EQ(true, tmp.persistence());
        for (size_type j=0; j < tmp.size(); ++j) {
            ASSERT_EQ(j&sdsl::bits::lo_set[exp_w_1], (size_type)tmp[j]);
        }
        tmp.persistence(false);
    }
}

TEST_F(IntVectorBufferTest, Swap)
{
    test_swap< sdsl::int_vector_buffer<> >(11, 2, vec_sizes);
    test_swap< sdsl::int_vector_buffer<1> >(1, 1, vec_sizes);
    test_swap< sdsl::int_vector_buffer<8> >(8, 8, vec_sizes);
    test_swap< sdsl::int_vector_buffer<16> >(16, 16, vec_sizes);
    test_swap< sdsl::int_vector_buffer<32> >(32, 32, vec_sizes);
    test_swap< sdsl::int_vector_buffer<64> >(64, 64, vec_sizes);
}

template<class t_T, class t_V>
void test_existing_file_persistence_close(size_type exp_w)
{
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;
    {
        // write file using int_vector_buffer
        t_T ivb(file_name, false, buffersize, 42, true);
        for (size_type i=0; i<1000000; ++i) {
            ivb[i] = rand();
        }
        ASSERT_EQ((size_type)1000000, ivb.size());
        ASSERT_EQ((size_type)exp_w, (size_type)ivb.width());
        ivb.close();
        // load int int_vector and int_vector_buffer and compare

        t_V iv;
        sdsl::load_from_file(iv, file_name);
        ASSERT_EQ((size_type)1000000, iv.size());
        t_T ivb2(file_name, true, buffersize, 5, false);
        ASSERT_EQ((size_type)1000000, ivb2.size());
        ASSERT_EQ((size_type)exp_w, (size_type)ivb2.width());
        for (size_type i=0; i<iv.size(); ++i) {
            ASSERT_EQ((size_type)ivb2[i], (size_type)iv[i]);
        }
        ivb2.close();
    }
}

//! Test the persistence in constructor, opening of existing file and close()
TEST_F(IntVectorBufferTest, ExistingFileAndPersistenceAndClose)
{
    test_existing_file_persistence_close< sdsl::int_vector_buffer<>, sdsl::int_vector<> >(42);
    test_existing_file_persistence_close< sdsl::int_vector_buffer<13>, sdsl::int_vector<13> >(13);
    test_existing_file_persistence_close< sdsl::int_vector_buffer<1>, sdsl::int_vector<1> >(1);
    test_existing_file_persistence_close< sdsl::int_vector_buffer<8>, sdsl::int_vector<8> >(8);
    test_existing_file_persistence_close< sdsl::int_vector_buffer<16>, sdsl::int_vector<16> >(16);
    test_existing_file_persistence_close< sdsl::int_vector_buffer<32>, sdsl::int_vector<32> >(32);
    test_existing_file_persistence_close< sdsl::int_vector_buffer<64>, sdsl::int_vector<64> >(64);
}

template<class t_T>
void test_add_and_sub(unsigned char w=1)
{
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 1024;
    t_T ivb(file_name, false, buffersize, w);
    for (size_type i=0; i < 10000; ++i) {
        value_type x = rand() & sdsl::bits::lo_set[ivb.width()];
        ivb[i] = x;
    }
    for (size_type i=0; i < ivb.size(); ++i) {
        value_type x = ivb[i];
        value_type y = rand() & sdsl::bits::lo_set[ivb.width()];
        ivb[i] += y;
        ASSERT_EQ((x+y)&sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ivb[i] -= y;
        ASSERT_EQ(x & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ivb[i]++;
        ASSERT_EQ((x+1)&sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ivb[i]--;
        ASSERT_EQ(x & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        ++ivb[i];
        ASSERT_EQ((x+1)&sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        --ivb[i];
        ASSERT_EQ(x & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
    }
}

TEST_F(IntVectorBufferTest, AddAndSub)
{
    for (unsigned char w=1; w <= 64; ++w) {
        test_add_and_sub< sdsl::int_vector_buffer<> >(w);
    }
    test_add_and_sub< sdsl::int_vector_buffer<1> >();
    test_add_and_sub< sdsl::int_vector_buffer<8> >();
    test_add_and_sub< sdsl::int_vector_buffer<16> >();
    test_add_and_sub< sdsl::int_vector_buffer<32> >();
    test_add_and_sub< sdsl::int_vector_buffer<64> >();
}

template<class t_T>
void test_random_access(unsigned char w=1)
{
    std::string file_name = "tmp/int_vector_buffer";
    size_type buffersize = 100;
    t_T ivb(file_name, false, buffersize, w);
    ASSERT_LT(ivb.buffersize()-buffersize, ivb.width());
    ASSERT_EQ((size_type)0, ivb.size());
    ivb[999] = 999 & sdsl::bits::lo_set[ivb.width()];
    ASSERT_EQ((size_type)999 & sdsl::bits::lo_set[ivb.width()], (size_type)ivb[999]);
    ASSERT_EQ((size_type)1000, ivb.size());
    for (size_type i=0; i < 999; ++i) {
        ASSERT_EQ((size_type)0, (size_type)ivb[i]);
    }
    for (size_type i=0; i < 1000; ++i) {
        value_type x = rand()%1000 & sdsl::bits::lo_set[ivb.width()];
        ivb[x] = x;
    }
    ASSERT_EQ((size_type)1000, ivb.size());
    for (size_type i=0; i < ivb.size(); ++i) {
        if (ivb[i]) {
            ASSERT_EQ(i&sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        }
    }
    buffersize = 50;
    ivb.buffersize(buffersize);
    ASSERT_EQ((size_type)1000, ivb.size());
    ASSERT_LT(ivb.buffersize()-buffersize, ivb.width());
    for (size_type i=0; i < 1000; ++i) {
        value_type x = rand()%1000 & sdsl::bits::lo_set[ivb.width()];
        ivb[x] = x;
    }
    for (size_type i=0; i < ivb.size(); ++i) {
        if (ivb[i]) {
            ASSERT_EQ(i&sdsl::bits::lo_set[ivb.width()], (size_type)ivb[i]);
        }
    }
}

TEST_F(IntVectorBufferTest, RandomAccess)
{
    for (unsigned char w=1; w <= 64; ++w) {
        test_random_access< sdsl::int_vector_buffer<> >(w);
    }
    test_random_access< sdsl::int_vector_buffer<1> >();
    test_random_access< sdsl::int_vector_buffer<8> >();
    test_random_access< sdsl::int_vector_buffer<16> >();
    test_random_access< sdsl::int_vector_buffer<32> >();
    test_random_access< sdsl::int_vector_buffer<64> >();
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
