#include "sdsl/nn_dict_dynamic.hpp"
#include "sdsl/nearest_neighbour_dictionary.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <string>
#include <random>

namespace
{

// The fixture for testing class nn_dict_dynamic.
class NNDictDynamicTest : public ::testing::Test
{
    protected:

        NNDictDynamicTest() {}

        virtual ~NNDictDynamicTest() {}

        virtual void SetUp() {}

        virtual void TearDown() {}
};

void compare_bv_and_nndd(const sdsl::bit_vector& bv, const sdsl::nn_dict_dynamic& nndd)
{
    sdsl::nearest_neighbour_dictionary<32> exp(bv);
    uint64_t first_one = exp.select(1);
    uint64_t last_one = exp.select(exp.rank(exp.size()));
    for (uint64_t i=0; i<first_one; ++i) {
        ASSERT_EQ(exp.size(), nndd.prev(i));
        ASSERT_EQ(exp.next(i), nndd.next(i));
    }
    for (uint64_t i=first_one; i<=last_one; ++i) {
        ASSERT_EQ(exp.prev(i), nndd.prev(i));
        ASSERT_EQ(exp.next(i), nndd.next(i));
    }
    for (uint64_t i=last_one+1; i<exp.size(); ++i) {
        ASSERT_EQ(exp.prev(i), nndd.prev(i));
        ASSERT_EQ(exp.size(), nndd.next(i));
    }
}

//! Test Constructors
TEST_F(NNDictDynamicTest, Constructors)
{
    static_assert(sdsl::util::is_regular<sdsl::nn_dict_dynamic>::value, "Type is not regular");
    uint64_t testsize = 100000;
    sdsl::bit_vector bv(testsize, 0);
    sdsl::nn_dict_dynamic nndd(testsize);

    // Fill nndd
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, testsize-1);
    auto dice = bind(distribution, rng);
    for (uint64_t i=0; i<testsize/4; ++i) {
        uint64_t value = dice();
        if (bv[value]) {
            bv[value] = 0;
            nndd[value] = 0;
        } else {
            bv[value] = 1;
            nndd[value] = 1;
        }
    }

    // Copy-constructor
    sdsl::nn_dict_dynamic nndd2(nndd);
    compare_bv_and_nndd(bv, nndd2);

    // Move-constructor
    sdsl::nn_dict_dynamic nndd3(std::move(nndd2));
    compare_bv_and_nndd(bv, nndd3);
    ASSERT_EQ((uint64_t)0, nndd2.size());

    // Copy-Assign
    sdsl::nn_dict_dynamic nndd4;
    nndd4 = nndd3;
    compare_bv_and_nndd(bv, nndd4);

    // Move-Assign
    sdsl::nn_dict_dynamic nndd5;
    nndd5 = std::move(nndd4);
    compare_bv_and_nndd(bv, nndd5);
    ASSERT_EQ((uint64_t)0, nndd4.size());
}

//! Test Operations next and prev
TEST_F(NNDictDynamicTest, NextAndPrev)
{
    uint64_t testsize = 100000;
    sdsl::bit_vector bv(testsize, 0);
    sdsl::nn_dict_dynamic nndd(testsize);
    for (uint64_t ones=1; ones<testsize; ones *= 2) {
        std::mt19937_64 rng(ones);
        std::uniform_int_distribution<uint64_t> distribution(0, testsize-1);
        auto dice = bind(distribution, rng);
        for (uint64_t i=0; i<ones; ++i) {
            uint64_t value = dice();
            if (bv[value]) {
                bv[value] = 0;
                nndd[value] = 0;
            } else {
                bv[value] = 1;
                nndd[value] = 1;
            }
        }
        bv[testsize/4] = 1;
        nndd[testsize/4] = 1;
        bv[3*testsize/4] = 1;
        nndd[3*testsize/4] = 1;
        compare_bv_and_nndd(bv, nndd);
    }
}

//! Test Serialize and Load
TEST_F(NNDictDynamicTest, SerializeAndLoad)
{
    std::string file_name = "tmp/nn_dict_dynamic";
    uint64_t testsize = 100000;
    sdsl::bit_vector bv(testsize, 0);
    {
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, testsize-1);
        auto dice = bind(distribution, rng);
        sdsl::nn_dict_dynamic nndd(testsize);
        for (uint64_t i=0; i<testsize/4; ++i) {
            uint64_t value = dice();
            if (bv[value]) {
                bv[value] = 0;
                nndd[value] = 0;
            } else {
                bv[value] = 1;
                nndd[value] = 1;
            }
        }
        sdsl::store_to_file(nndd, file_name);
    }
    {
        sdsl::nn_dict_dynamic nndd(0);
        sdsl::load_from_file(nndd, file_name);
        compare_bv_and_nndd(bv, nndd);
    }
    sdsl::remove(file_name);
}



}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
