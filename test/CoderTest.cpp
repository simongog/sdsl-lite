#include "sdsl/int_vector.hpp"
#include "sdsl/coder.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <random>

namespace
{

using namespace sdsl;

// The fixture for testing class int_vector.
template<class T>
class CoderTest : public ::testing::Test
{
    protected:

        CoderTest() {
            m_data = sdsl::int_vector<>(10000000);
            util::set_random_bits(m_data);
            for (size_t i=0; i<m_data.size()/3; ++i)
                m_data[i] = i;

            std::mt19937_64 rng;
            std::uniform_int_distribution<uint64_t> distribution(0, 6);
            auto dice = bind(distribution, rng);
            for (size_t i=m_data.size()/3; i<2*m_data.size()/3; ++i)
                m_data[i] = dice();
        }

        virtual ~CoderTest() { }

        virtual void SetUp() { }

        virtual void TearDown() { }
        sdsl::int_vector<> m_data;
};

using testing::Types;
typedef Types<
      coder::elias_delta,
      coder::elias_gamma,
      coder::fibonacci,
      coder::comma<>,
      coder::comma<4>,
      coder::comma<8>,
      coder::comma<16>
      >
      Implementations;

TYPED_TEST_CASE(CoderTest, Implementations);

TYPED_TEST(CoderTest, SinlgeEncodeDecode)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    uint8_t offset = 0;
    uint64_t buf[8] = {0};
    uint64_t* pb = buf;
    for (size_t i=0; i<this->m_data.size(); ++i) {
        TypeParam::encode(this->m_data[i], pb, offset);
        pb = buf;
        offset = 0;
        uint64_t x = TypeParam::template decode<false,false,uint64_t*>(buf, 0, 1);
        ASSERT_EQ(this->m_data[i], x);
    }
}

TYPED_TEST(CoderTest, AllEncodeDecode)
{
    int_vector<> tmp,data;
    TypeParam::encode(this->m_data, tmp);
    TypeParam::decode(tmp, data);
    ASSERT_EQ(this->m_data.size(), data.size());
    for (size_t i=0; i<this->m_data.size(); ++i) {
        ASSERT_EQ(this->m_data[i], data[i]);
    }
}

TYPED_TEST(CoderTest, DecodePrefixSum)
{
    int_vector<> tmp;
    TypeParam::encode(this->m_data, tmp);
    uint64_t start = 0;
    const uint64_t sample = 32;
    for (size_t i=0; i<this->m_data.size(); ++i) {
        uint64_t sum = 0;
        uint64_t lstart = start;
        if (i%sample==0) {
            for (size_t j=1; j<=sample and i+j-1<this->m_data.size(); ++j) {
                sum += this->m_data[i+j-1];
                ASSERT_EQ(sum, TypeParam::decode_prefix_sum(tmp.data(), lstart, j));
                uint64_t x = TypeParam::template decode<true,false,uint64_t*>(tmp.data(), lstart, j);
                ASSERT_EQ(sum, x);
            }
        }
        start += TypeParam::encoding_length(this->m_data[i]);
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
