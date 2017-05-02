#include <algorithm>
#include <chrono>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "sdsl/dac_vector.hpp"
#include "sdsl/util.hpp"

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::int_vector<>::value_type value_type;

void run_test(const std::vector<value_type>& vec) {
    {
        sdsl::dac_vector<> v(vec);
        std::stringstream ss;
        v.serialize(ss);
        std::cout << "  old = " << ss.str().size() << std::endl;
        sdsl::dac_vector<> w;
        w.load(ss);

        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], w[i]);
    }
    size_t last_levels = -1;
    for (int max_levels = 1; max_levels <= 32; max_levels *= 2) {
        sdsl::dac_vector_dp<> v(vec, max_levels);
        if (v.levels() == last_levels)
            break;
        last_levels = v.levels();

        std::stringstream ss;
        v.serialize(ss);
        std::cout << "  new (plain) with " << v.levels() << " levels = "
            << ss.str().size() << std::endl;
        sdsl::dac_vector_dp<> w;
        w.load(ss);
        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], w[i]);

        // test move
        sdsl::dac_vector_dp<> z = std::move(w);
        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], z[i]);

        // test copy
        w = z;
        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], w[i]);
    }

    for (int max_levels = 1; max_levels <= 32; max_levels *= 2) {
        sdsl::dac_vector_dp<sdsl::rrr_vector<>> v(vec, max_levels);
        if (v.levels() == last_levels)
            break;
        last_levels = v.levels();

        std::stringstream ss;
        v.serialize(ss);
        std::cout << "  new (rrr) with " << v.levels() << " levels = "
            << ss.str().size() << std::endl;
        sdsl::dac_vector_dp<sdsl::rrr_vector<>> w;
        w.load(ss);
        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], w[i]);

        // test move
        sdsl::dac_vector_dp<> z = std::move(w);
        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], z[i]);

        // test copy
        w = z;
        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], w[i]);
    }
}

TEST(DacVectorTest, SmokeTest) {
    std::vector<value_type> vec { 1, 3, 5, 7, 100, 2000000000 };
    sdsl::dac_vector_dp<> v(vec, 3);
    ASSERT_EQ(3ul, v.levels());
    for (size_t i = 0; i < vec.size(); ++i)
        ASSERT_EQ(vec[i], v[i]);
}

TEST(DacVectorTest, SmokeTestLarge) {
    std::vector<value_type> vec, add { 1, 3, 5, 7, 100, 2000000000 };
    for (int i = 0; i < 50; ++i)
        std::copy(add.begin(), add.end(), std::back_inserter(vec));
    run_test(vec);
}

TEST(DacVectorTest, LargeRandom) {
    for (int i = 0; i < 10; ++i) {
        std::cout << "random uniform test " << i << std::endl;
        std::vector<value_type> vec;
        std::mt19937_64 rng;
        {
            std::uniform_int_distribution<uint64_t> distribution(0, 100000000);
            auto dice = bind(distribution, rng);
            for (size_t i=0; i < 100000; ++i)
                vec.push_back(dice());
        }
        run_test(vec);
    }

    for (int i = 0; i < 10; ++i) {
        std::cout << "random exponential test " << i << std::endl;
        std::vector<value_type> vec;
        std::mt19937_64 rng;
        {
            std::exponential_distribution<double> dist(3.5);
            auto dice = bind(dist, rng);
            for (size_t i = 0; i < 100000; ++i)
                vec.push_back(uint64_t(dice() * 1000000));
        }
        run_test(vec);
    }
}

template <typename F>
uint64_t measure(F f) {
    using timer = std::chrono::high_resolution_clock;
    auto start = timer::now();
    f();
    return std::chrono::duration_cast<std::chrono::microseconds>(
            timer::now() - start).count();
}

TEST(DacVectorTest, LargeRandomBench) {
    for (size_t sz : {10000, 100000, 1000000, 10000000}) {
        size_t queries = 100000;
        std::vector<value_type> vec;
        std::vector<size_t> query_indexes;

        std::mt19937_64 rng;
        {
            std::exponential_distribution<double> dist(3.5);
            auto dice = bind(dist, rng);
            for (size_t i=0; i < sz; ++i)
                vec.push_back(uint64_t(dice() * 1000000));
        }

        {
            std::uniform_int_distribution<uint64_t> distribution(0, sz - 1);
            auto dice = bind(distribution, rng);
            for (size_t i=0; i < queries; ++i)
                query_indexes.push_back(dice());
        }

        sdsl::dac_vector<> v1(vec);
        sdsl::dac_vector_dp<> v2(vec);
        sdsl::dac_vector_dp<sdsl::rrr_vector<>> v3(vec);

        std::cout << "benchmark size " << sz << std::endl;
        std::cout << "  dac_vector levels    " << v1.levels() << std::endl;
        std::cout << "  dac_vector_dp levels " << v2.levels() << std::endl;

        size_t s1=0;
        std::cout << "  time dac_vector                  " << measure([&]() {
            for (size_t query : query_indexes)
                s1+=v1[query];
        }) << std::endl;

        size_t s2=0;
        std::cout << "  time dac_vector_dp               " << measure([&]() {
            for (size_t query : query_indexes)
                s2+=v2[query];
        }) << std::endl;

        size_t s3=0;
        std::cout << "  time dac_vector_dp<rrr_vector<>> " << measure([&]() {
            for (size_t query : query_indexes)
                s3+=v3[query];
        }) << std::endl;

        ASSERT_EQ(s1, s2);
        ASSERT_EQ(s1, s3);
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
