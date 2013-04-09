#include "sdsl/bit_vectors.hpp" // for rrr_vector
#include "gtest/gtest.h"
#include <string>

using namespace sdsl;
using namespace std;

string test_file;

namespace
{

template<class T>
class BitVectorTest : public ::testing::Test { };

using testing::Types;

typedef Types<
bit_vector,
bit_vector_il<256>,
bit_vector_il<512>,
bit_vector_il<1024>,
rrr_vector<64>,
rrr_vector<256>,
rrr_vector<129>,
rrr_vector<192>,
rrr_vector<255>,
rrr_vector<15>,
rrr_vector<31>,
rrr_vector<63>,
rrr_vector<83>,
rrr_vector<127>,
rrr_vector<128>,
sd_vector<>,
sd_vector<rrr_vector<63> >
> Implementations;

TYPED_TEST_CASE(BitVectorTest, Implementations);

//! Test operator[]
TYPED_TEST(BitVectorTest, Access)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    for (uint64_t j=0; j < bv.size(); ++j) {
        ASSERT_EQ((bool)(bv[j]), (bool)(c_bv[j]));
    }
}

TYPED_TEST(BitVectorTest, Swap)
{
    bit_vector bv;
    ASSERT_TRUE(load_from_file(bv, test_file));
    TypeParam c_bv(bv);
    ASSERT_EQ(bv.size(), c_bv.size());
    TypeParam bv_empty;
    ASSERT_EQ((uint64_t)0, bv_empty.size());
    bv_empty.swap(c_bv);
    ASSERT_EQ((uint64_t)0, c_bv.size());
    ASSERT_EQ(bv.size(), bv_empty.size());
    for (uint64_t j=0; j < bv.size(); ++j) {
        ASSERT_EQ((bool)(bv[j]), (bool)(bv_empty[j]));
    }
}

}// end namespace

int main(int argc, char* argv[])
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

