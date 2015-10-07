#include "sdsl/sd_vector.hpp"
#include "gtest/gtest.h"

using namespace sdsl;
using namespace std;

namespace
{

const size_t BV_SIZE = 1000000;

TEST(sd_vector_test, iterator_constructor)
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
    sd_vector<> sdv(pos.begin(),pos.end());
    for (size_t i=0; i < bv.size(); ++i) {
        ASSERT_EQ((bool)sdv[i],(bool)bv[i]);
    }
}

TEST(sd_vector_test, builder_constructor)
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
    sd_vector<> sdv(builder);
    for (size_t i=0; i < bv.size(); ++i) {
        ASSERT_EQ((bool)sdv[i],(bool)bv[i]);
    }
}

} // end namespace

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

