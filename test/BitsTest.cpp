#include "sdsl/int_vector.hpp"
#include "sdsl/bits.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>

namespace
{

SDSL_UNUSED uint32_t cnt_naive(uint64_t x)
{
    uint32_t res = 0;
    for (int i=0; i<64; ++i) {
        res += x&1;
        x>>=1;
    }
    return res;
}

SDSL_UNUSED inline uint32_t cnt10_naive(uint64_t x, uint64_t& c)
{
    uint32_t sum = 0, lastbit = c;
    for (uint32_t i=0; i<64; ++i) {
        if ((x&1) == 0 and lastbit == 1) {
            ++sum;
        }
        lastbit = (x&1);
        x >>= 1;
    }
    c = lastbit;
    return sum;
}

SDSL_UNUSED uint32_t cnt11_naive(uint64_t x)
{
    uint32_t res = 0;
    for (int i=0; i<64; ++i) {
        if ((x&3) == 3) {
            ++res;
            ++i;
            x>>=1;
        }
        x>>=1;
    }
    return res;
}

SDSL_UNUSED inline uint32_t hi_naive(uint64_t x)
{
    uint32_t res=63;
    while (x) {
        if (x&0x8000000000000000ULL) {
            return res;
        }
        --res;
        x<<=1;
    }
    return 0;
}

inline uint32_t lo_naive(uint64_t x)
{
    if (x&1)
        return 0;
    x>>=1;
    for (int i=1; i<64; ++i)
        if (x&1)
            return i;
        else
            x>>=1;
    return 63;
}

SDSL_UNUSED inline uint32_t sel_naive(uint64_t x, uint32_t i)
{
    uint32_t pos = 0;
    while (x) {
        i -= x&1;
        x>>=1;
        if (!i) break;
        ++pos;
    }
    return pos;
}

SDSL_UNUSED uint32_t sel11_naive(uint64_t x, uint32_t i)
{
    for (uint32_t j=0; j<63; ++j) {
        if ((x&3)==3) {
            i--;
            if (!i) return j+1;
            x>>=1;
            ++j;
        }
        x>>=1;
    }
    return 63;
}

SDSL_UNUSED inline uint64_t rev_naive(uint64_t x)
{
    uint64_t y = 0;
    for (size_t i=0; i<64; i++) {
        if (x&(1ULL << i)) {
            y |= (1ULL << (63-i));
        }
    }
    return y;
}


// The fixture for testing class int_vector.
class BitsTest : public ::testing::Test
{
    protected:

        BitsTest() {
            m_data = sdsl::int_vector<64>(1000000);
            sdsl::util::set_random_bits(m_data);
        }

        virtual ~BitsTest() { }

        virtual void SetUp() { }

        virtual void TearDown() { }
        sdsl::int_vector<64> m_data;
};

//! Test the default constructor
TEST_F(BitsTest, cnt)
{
    for (uint64_t i=0; i<64; ++i) {
        ASSERT_EQ((uint64_t)1, sdsl::bits::cnt(1ULL<<i));
    }
}

//! Test the parametrized constructor
TEST_F(BitsTest, sel)
{
    for (uint64_t i=0; i<64; ++i) {
        ASSERT_EQ(i, sdsl::bits::sel(1ULL<<i, 1));
    }
    for (uint64_t i=0; i < this->m_data.size(); ++i) {
        uint64_t x = this->m_data[i];
        for (uint64_t j=0, ones=0; j<64; ++j) {
            if ((x >> j)&1) {
                ++ones;
                ASSERT_EQ(j, sdsl::bits::sel(x, ones));
            }
        }
    }
}

TEST_F(BitsTest, hi)
{
    for (uint64_t i=0; i<64; ++i) {
        uint64_t x = 1ULL<<i;
        ASSERT_EQ(hi_naive(x), sdsl::bits::hi(x));
    }
}


TEST_F(BitsTest, lo)
{
    for (uint64_t i=0; i<64; ++i) {
        uint64_t x = 1ULL<<i;
        ASSERT_EQ(lo_naive(x), sdsl::bits::lo(x));
    }
}

TEST_F(BitsTest, rev)
{
    for (uint64_t i=0; i < this->m_data.size(); ++i) {
        uint64_t x = this->m_data[i];
        uint64_t rx = sdsl::bits::rev(x);
        ASSERT_EQ(rev_naive(x),rx);
        ASSERT_EQ(x,sdsl::bits::rev(rx));
    }
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
