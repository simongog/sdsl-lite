#include <algorithm>
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
        sdsl::dac_vector_dp v(vec);
        std::stringstream ss;
        v.serialize(ss);
        std::cout << "new = " << ss.str().size() << std::endl;
        sdsl::dac_vector_dp w;
        w.load(ss);

        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], w[i]);
    }
    {
        sdsl::dac_vector<> v(vec);
        std::stringstream ss;
        v.serialize(ss);
        std::cout << "old = " << ss.str().size() << std::endl;
        sdsl::dac_vector<> w;
        w.load(ss);

        for (size_t i = 0; i < vec.size(); ++i)
            ASSERT_EQ(vec[i], w[i]);
    }
}

TEST(DacVectorTest, SmokeTest) {
    std::vector<value_type> vec, add { 1, 3, 5, 7, 100, 2000000000 };
    for (int i = 0; i < 50; ++i)
        std::copy(add.begin(), add.end(), std::back_inserter(vec));
    run_test(vec);
}

TEST(DacVectorTest, LargeRandom) {
    for (int i = 0; i < 10; ++i) {
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

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
