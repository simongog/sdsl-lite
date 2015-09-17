#include "sdsl/bit_vectors.hpp"
#include "sdsl/rank_support.hpp"
#include "gtest/gtest.h"
#include <string>

using namespace sdsl;
using namespace std;

string test_file;

namespace
{

template<class T>
class rank_support_test : public ::testing::Test { };

using testing::Types;

typedef Types<rank_support_il<1, 256>,
        rank_support_il<1, 512>,
        rank_support_il<1, 1024>,
        rank_support_rrr<>,
        rank_support_v<>,
        rank_support_v5<>,
        rank_support_rrr<1, 64>,
        rank_support_rrr<1, 192>,
        rank_support_rrr<1, 256>,
        rank_support_rrr<1, 255>,
        rank_support_rrr<1, 15>,
        rank_support_rrr<1, 31>,
        rank_support_rrr<1, 63>,
        rank_support_rrr<1, 83>,
        rank_support_rrr<1, 127>,
        rank_support_rrr<1, 128>,
        rank_support_rrr<1, 129>,
        rank_support_sd<1>,
        rank_support_il<0, 256>,
        rank_support_il<0, 512>,
        rank_support_il<0, 1024>,
        rank_support_rrr<0>,
        rank_support_v<0>,
        rank_support_v5<0>,
        rank_support_rrr<0, 64>,
        rank_support_rrr<0, 192>,
        rank_support_rrr<0, 256>,
        rank_support_rrr<0, 255>,
        rank_support_rrr<0, 15>,
        rank_support_rrr<0, 30>,
        rank_support_rrr<0, 63>,
        rank_support_rrr<0, 83>,
        rank_support_rrr<0, 127>,
        rank_support_rrr<0, 128>,
        rank_support_rrr<0, 129>,
        rank_support_sd<0>,
        rank_support_hyb<1>,
        rank_support_hyb<0>,
        rank_support_v<10,2>,
        rank_support_v<01,2>,
        rank_support_v<00,2>,
        rank_support_v<11,2>,
        rank_support_v5<10,2>,
        rank_support_v5<01,2>,
        rank_support_v5<00,2>,
        rank_support_v5<11,2>
        > Implementations;

TYPED_TEST_CASE(rank_support_test, Implementations);

//! Test the rank method
TYPED_TEST(rank_support_test, rank_method)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bvec;
    ASSERT_TRUE(load_from_file(bvec, test_file));
    typename TypeParam::bit_vector_type bv(bvec);
    TypeParam rs(&bv);
    uint64_t rank=0;
    for (uint64_t j=0; j < bvec.size(); ++j) {
        ASSERT_EQ(rank, rs.rank(j));
        bool found = (j >= TypeParam::bit_pat_len-1);
        for (uint8_t k=0; found and k < TypeParam::bit_pat_len; ++k) {
            found &= bvec[j-k] == ((TypeParam::bit_pat>>k)&1);
        }
        rank += found;
    }
    EXPECT_EQ(rank, rs.rank(bvec.size()));
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " FILE " << endl;
        cout << "  Reads a bitvector from FILE and executes tests." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    return RUN_ALL_TESTS();
}

