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

class SDVectorTest : public ::testing::Test { };

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

TEST_F(SDVectorTest, getint)
{
    ASSERT_EQ(1+1,2);
    bit_vector bv(10000, 0);
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, 9);
    auto dice = bind(distribution, rng);
    for (size_t i=0; i < bv.size(); ++i) {
        if (0 == dice())
            bv[i] = 1;
    }
    for (size_t i=0; i<= 1000; ++i)
        bv[i] = 0;

    sd_vector<> sdb(bv);
    for (size_t len=1; len<=64; ++len) {
        for (size_t i=0; i+len <= bv.size(); ++i) {
            ASSERT_EQ(bv.get_int(i,len), sdb.get_int(i,len))
                    << "i="<<i<<" len="<<len<<endl;
        }
    }
}

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
    TypeParam mo_bv = TypeParam(bv);
    ASSERT_EQ(bv.size(), mo_bv.size());
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
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " FILE " << endl;
        cout << "  Reads a bitvector from FILE and executes tests." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    return RUN_ALL_TESTS();
}

