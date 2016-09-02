#include "sdsl/sd_vector.hpp"
#include "sdsl/bit_vectors.hpp"
#include "gtest/gtest.h"

using namespace sdsl;
using namespace std;

namespace
{

const size_t BV_SIZE = 1000000;

template<class T>
class sd_vector_test : public ::testing::Test { };

using testing::Types;

typedef Types<
sd_vector<>,
          sd_vector<rrr_vector<63>>
          > Implementations;

TYPED_TEST_CASE(sd_vector_test, Implementations);

TYPED_TEST(sd_vector_test, iterator_constructor)
{
    std::vector<uint64_t> pos;
    bit_vector bv(BV_SIZE);
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, 9);
    auto dice = bind(distribution, rng);
    for (size_t i=0; i < bv.size(); ++i) {
        if (0 == dice()) {
            pos.emplace_back(i);
            bv[i] = 1;
        }
    }
    TypeParam sdv(pos.begin(),pos.end());
    for (size_t i=0; i < bv.size(); ++i) {
        ASSERT_EQ((bool)sdv[i],(bool)bv[i]);
    }
}

TYPED_TEST(sd_vector_test, builder_constructor)
{
    std::vector<uint64_t> pos;
    bit_vector bv(BV_SIZE);
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, 9);
    auto dice = bind(distribution, rng);
    size_t ones = 0;
    for (size_t i=0; i < bv.size(); ++i) {
        if (0 == dice()) {
            pos.emplace_back(i);
            bv[i] = 1;
            ones++;
        }
    }
    sd_vector_builder builder(BV_SIZE, ones);
    for (auto i : pos) {
        builder.set(i);
    }
    TypeParam sdv(builder);
    for (size_t i=0; i < bv.size(); ++i) {
        ASSERT_EQ((bool)sdv[i],(bool)bv[i]);
    }
}

TYPED_TEST(sd_vector_test, builder_empty_constructor)
{
    sd_vector_builder builder(BV_SIZE, 0UL);
    TypeParam sdv(builder);
    for (size_t i=0; i < BV_SIZE; ++i) {
        ASSERT_FALSE((bool)sdv[i]);
    }
}

} // end namespace

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

