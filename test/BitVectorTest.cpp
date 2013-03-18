#include "sdsl/bit_vectors.hpp" // for rrr_vector
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::bit_vector bit_vector;

template<class T>
class BitVectorTest : public ::testing::Test
{
    protected:

        static const size_t n = 40;

        BitVectorTest() {
            // You can do set-up work for each test here
        }

        virtual ~BitVectorTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            srand(13);
            // crafted small examples
            bs[0] = bit_vector(32,0);
            bs[0][1] = bs[0][4] = bs[0][7] = bs[0][18] =
                                                 bs[0][24] = bs[0][26] = bs[0][30] = bs[0][31] = 1;
            bs[1] = bit_vector(1,0);
            bs[2] = bit_vector(1000000,0);
            bs[3] = bit_vector(1000000,1);
            bs[4] = bit_vector(0);
            for (size_type i=5; i<14; ++i) {
                bs[i] = bit_vector(i, i%2);
            }
            // dense vectors of different sizes
            for (size_type i=14; i<n-4; ++i) {
                bs[i] = bit_vector(rand() % (8<< (i-14)));
                for (size_type j=0; j < bs[i].size(); ++j) {
                    if (rand() % 2)
                        bs[i][j] = 1;
                }
            }
            size_type last_size = 1000000;
            bs[n-4] = bit_vector(last_size, 1);
            bs[n-3] = bit_vector(last_size, 0);
            // populate vectors with some other bits
            for (size_type i=0; i<last_size/1000; ++i) {
                size_type x = rand()%last_size;
                bs[n-4][x] = 0;
                bs[n-3][x] = 1;
            }
            bs[n-2] = bit_vector(last_size, 1);
            bs[n-1] = bit_vector(last_size, 0);
            // populate vectors with some blocks of other bits
            for (size_type i=0; i<last_size/1000; ++i) {
                size_type x = rand()%last_size;
                size_type len = rand()%1000;
                for (size_type j=x; j<x+len and j<last_size; ++j) {
                    bs[n-2][j] = 0;
                    bs[n-1][j] = 1;
                }
            }
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test case for Foo.
        bit_vector bs[n];
};

using testing::Types;

typedef Types<
sdsl::bit_vector,
     sdsl::bit_vector_il<256>,
     sdsl::bit_vector_il<512>,
     sdsl::bit_vector_il<1024>,
     sdsl::rrr_vector<64>,
     sdsl::rrr_vector<256>,
     sdsl::rrr_vector<129>,
     sdsl::rrr_vector<192>,
     sdsl::rrr_vector<255>,
     sdsl::rrr_vector<15>,
     sdsl::rrr_vector<31>,
     sdsl::rrr_vector<63>,
     sdsl::rrr_vector<83>,
     sdsl::rrr_vector<127>,
     sdsl::rrr_vector<128>,
     sdsl::sd_vector<>,
     sdsl::sd_vector<sdsl::rrr_vector<63> >
     > Implementations;

TYPED_TEST_CASE(BitVectorTest, Implementations);

//! Test operator[]
TYPED_TEST(BitVectorTest, Access)
{
    for (size_type i=0; i<this->n; ++i) {
        TypeParam copied_bs(this->bs[i]);
        ASSERT_EQ((this->bs[i]).size(), copied_bs.size());
        for (size_type j=0; j < (this->bs[i]).size(); ++j) {
            ASSERT_EQ((bool)(this->bs[i][j]), (bool)(copied_bs[j])) << " at index "<<j<<" of vector "<<i<<" of length "<<(this->bs[i]).size();
        }
    }
}

TYPED_TEST(BitVectorTest, Swap)
{
    for (size_type i=0; i<this->n; ++i) {
        TypeParam copied_bs(this->bs[i]);
        ASSERT_EQ(this->bs[i].size(), copied_bs.size());

        TypeParam bs_empty;
        ASSERT_EQ((size_type)0, bs_empty.size());

        bs_empty.swap(copied_bs);
        ASSERT_EQ((size_type)0, copied_bs.size());
        ASSERT_EQ(this->bs[i].size(), bs_empty.size());
        for (size_type j=0; j < (this->bs[i]).size(); ++j) {
            ASSERT_EQ((bool)(this->bs[i][j]), (bool)(bs_empty[j])) << " at index "<<j<<" of vector "<<i<<" of length "<<(this->bs[i]).size();
        }
    }
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

