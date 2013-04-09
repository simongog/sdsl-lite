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

typedef Types<sdsl::select_support_mcl<0>,
        sdsl::rrr_select_support<0, 256>,
        sdsl::rrr_select_support<0>,
        sdsl::rrr_select_support<0, 15>,
        sdsl::rrr_select_support<0, 31>,
        sdsl::rrr_select_support<0, 63>,
        sdsl::rrr_select_support<0, 127>,
        sdsl::select_support_interleaved<0, 256>,
        sdsl::select_support_interleaved<0, 512>,
        sdsl::select_support_interleaved<0, 1024>
        > Implementations;

TYPED_TEST_CASE(SelectSupportTest, Implementations);

//! Test the rank method
TYPED_TEST(SelectSupportTest, SelectMethod)
{
    bit_vector bvec;
    ASSERT_TRUE(load_from_file(bvec, test_file));
    typename TypeParam::bit_vector_type bv(bvec);
    TypeParam ss(&bv);
    for (size_type j=0, select=0; j < (this->bs[i]).size(); ++j) {
        if (!bvec[j]) {
            ++select;
            ASSERT_EQ(ss.select(select), j);
        }
    }
}

}// end namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " FILE " << endl;
        cout << "  Reads a bitvector from FILE and executes tests." << endl;
        return 1;
    }
    test_file = argv[1];
    return RUN_ALL_TESTS();
}
