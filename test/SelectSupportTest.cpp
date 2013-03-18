
#include "sdsl/int_vector.hpp"
#include "sdsl/select_support_mcl.hpp" // for select_support_mcl
#include "sdsl/bit_vector_il.hpp" // for rank_support_interleaved
#include "sdsl/rrr_vector.hpp" // for select_support_rrr
#include "sdsl/sd_vector.hpp" // for select_support_sd
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::bit_vector bit_vector;

template<class T>
class SelectSupportTest : public ::testing::Test
{
    protected:

        static const size_t n = 41;

        SelectSupportTest() {
            // You can do set-up work for each test here
        }

        virtual ~SelectSupportTest() {
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
            for (size_type i=14; i<n-5; ++i) {
                bs[i] = bit_vector(rand() % (8<< (i-14)));
                for (size_type j=0; j < bs[i].size(); ++j) {
                    if (rand() % 2)
                        bs[i][j] = 1;
                }
            }
            size_type last_size = 1000000;
            bs[n-5] = bit_vector(last_size, 1);
            bs[n-4] = bit_vector(last_size, 0);
            // populate vectors with some other bits
            for (size_type i=0; i<last_size/1000; ++i) {
                size_type x = rand()%last_size;
                bs[n-5][x] = 0;
                bs[n-4][x] = 1;
            }
            bs[n-3] = bit_vector(last_size, 1);
            bs[n-2] = bit_vector(last_size, 0);
            // populate vectors with some blocks of other bits
            for (size_type i=0; i<last_size/1000; ++i) {
                size_type x = rand()%last_size;
                size_type len = rand()%1000;
                for (size_type j=x; j<x+len and j<last_size; ++j) {
                    bs[n-3][j] = 0;
                    bs[n-2][j] = 1;
                }
            }
            {
                // Matthias Petri's test
                srand(4711);
                size_type size = 1000000;
                size_type ones = 4030+(rand()%80);
                sdsl::util::assign(bs[n-1], bit_vector(size));
                for (size_type i=0; i<ones; i++) {
                    size_type rnd = rand();
                    while (bs[n-1][rnd%size]==1) {
                        rnd = rand();
                    }
                    bs[n-1][rnd%size] = 1;
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

typedef Types<sdsl::select_support_mcl<>,
        sdsl::select_support_rrr<1, 256>,
        sdsl::select_support_rrr<1, 129>,
        sdsl::select_support_rrr<1, 192>,
        sdsl::select_support_rrr<1, 255>,
        sdsl::select_support_rrr<1, 15>,
        sdsl::select_support_rrr<1, 31>,
        sdsl::select_support_rrr<1, 63>,
        sdsl::select_support_rrr<1, 127>,
        sdsl::select_support_rrr<1, 128>,
        sdsl::select_support_sd<>,
        sdsl::select_support_il<1, 256>,
        sdsl::select_support_il<1, 512>,
        sdsl::select_support_il<1, 1024>
        > Implementations;

TYPED_TEST_CASE(SelectSupportTest, Implementations);

//! Test the select method
TYPED_TEST(SelectSupportTest, SelectMethod)
{
    for (size_type i=0; i<this->n; ++i) {
        typename TypeParam::bit_vector_type bv(this->bs[i]);
        TypeParam ss(&bv);
        for (size_type j=0, select=0; j < (this->bs[i]).size(); ++j) {
            if (this->bs[i][j]) {
                ++select;
//				EXPECT_EQ(ss.select(select), j) << " at query "<<select<<" of vector "<<i<<" of length "<<(this->bs[i]).size();
                ASSERT_EQ(ss.select(select), j) << " at query "<<select<<" of vector "<<i<<" of length "<<(this->bs[i]).size();
            }
        }
    }
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

