#include "sdsl/bit_vectors.hpp"
#include "sdsl/select_support.hpp"
#include "gtest/gtest.h"
#include <string>

using namespace sdsl;
using namespace std;

string test_file;

namespace
{

template<class T>
class SelectSupportTest : public ::testing::Test { };

using testing::Types;

typedef Types<select_support_mcl<>,
        select_support_rrr<1, 256>,
        select_support_rrr<1, 129>,
        select_support_rrr<1, 192>,
        select_support_rrr<1, 255>,
        select_support_rrr<1, 15>,
        select_support_rrr<1, 31>,
        select_support_rrr<1, 63>,
        select_support_rrr<1, 127>,
        select_support_rrr<1, 128>,
        select_support_sd<1>,
        select_support_sd<0>,
        select_0_support_sd<>,
        select_support_il<1, 256>,
        select_support_il<1, 512>,
        select_support_il<1, 1024>,
        select_support_mcl<0>,
        select_support_rrr<0, 256>,
        select_support_rrr<0>,
        select_support_rrr<0, 15>,
        select_support_rrr<0, 31>,
        select_support_rrr<0, 63>,
        select_support_rrr<0, 127>,
        select_support_il<0, 256>,
        select_support_il<0, 512>,
        select_support_il<0, 1024>
        > Implementations;

TYPED_TEST_CASE(SelectSupportTest, Implementations);

//! Test the select method
TYPED_TEST(SelectSupportTest, SelectMethod)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    bit_vector bvec;
    ASSERT_TRUE(load_from_file(bvec, test_file));
    typename TypeParam::bit_vector_type bv(bvec);
    TypeParam ss(&bv);
    for (uint64_t j=0, select=0; j < bvec.size(); ++j) {
        if (bvec[j] == TypeParam::bit_pat) {
            ++select;
            ASSERT_EQ(j, ss.select(select));
        }
    }
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
