#include "sdsl/wavelet_trees.hpp"
#include "sdsl/rrr_vector.hpp"
#include "sdsl/bit_vector_il.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <algorithm> // for std::min

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;

string test_file;
string temp_file;

template<class T>
class WtByteTest : public ::testing::Test { };

using testing::Types;

typedef Types<
wt<unsigned char*, rrr_vector<63> >,
   wt<unsigned char*, bit_vector_il<>  >,
   wt<unsigned char*, bit_vector>,
   wt_huff<bit_vector_il<> >,
   wt_huff<bit_vector, rank_support_v<> >,
   wt_huff<bit_vector, rank_support_v5<> >,
   wt_huff<rrr_vector<63> >,
   wt_rlmn<>,
   wt_rlmn<bit_vector>,
   wt_rlg<>,
   wt_rlg8<>
   > Implementations;

TYPED_TEST_CASE(WtByteTest, Implementations);

TYPED_TEST(WtByteTest, CreateAndStoreTest)
{
    TypeParam wt;
    construct(wt, test_file, 1);
    bool success = store_to_file(wt, temp_file);
    ASSERT_EQ(true, success);
}

//! Test access methods
TYPED_TEST(WtByteTest, Sigma)
{
    TypeParam wt;
    ASSERT_EQ(true, load_from_file(wt, temp_file));
    int_vector<8> text;
    ASSERT_EQ(true, load_vector_from_file(text, test_file, 1));
    ASSERT_EQ(text.size(), wt.size());
    bit_vector occur(256, 0);
    uint16_t sigma = 0;
    for (size_type j=0; j<text.size(); ++j) {
        if (!occur[(unsigned char)text[j]]) {
            occur[(unsigned char)text[j]] = 1;
            ++sigma;
        }
    }
    ASSERT_EQ(sigma, wt.sigma);
}

//! Test access methods
TYPED_TEST(WtByteTest, Access)
{
    TypeParam wt;
    ASSERT_EQ(true, load_from_file(wt, temp_file));
    int_vector<8> text;
    ASSERT_EQ(true, load_vector_from_file(text, test_file, 1));
    ASSERT_EQ(text.size(), wt.size());
    for (size_type j=0; j<text.size(); ++j) {
        ASSERT_EQ((typename TypeParam::value_type)text[j], wt[j])<<" j="<<j;
    }
}

template<class tWt>
void test_rank(const tWt& wt, const int_vector<8>& text, size_type n)
{
    vector<size_type> cnt(256, 0);
    ASSERT_EQ(n, wt.size());
    for (size_type j=0; j < wt.size(); ++j) {
        cnt[text[j]]++;
        ASSERT_EQ(cnt[text[j]], wt.rank(j+1, text[j]))<< " j = "<<j<<" text[j]"<<text[j];
    }
    // Do random queries for all characters that do not occur in the string
    for (size_type j=0; j<cnt.size(); ++j) {
        if (cnt[j] == 0) {
            for (size_type k=0; k<1000; ++k) {
                size_type pos = rand()%(wt.size()+1);
                ASSERT_EQ((size_type)0, wt.rank(pos, (unsigned char)j))<<" pos="<<pos;
            }
        }
    }
    // Test rank(size(), c) for each character c
    for (size_type c=0; c < 256; ++c) {
        ASSERT_EQ(cnt[c], wt.rank(wt.size(), (unsigned char)c))<<" c="<<c;
    }
}

//! Test rank methods
TYPED_TEST(WtByteTest, Rank)
{
    TypeParam wt;
    ASSERT_EQ(true, load_from_file(wt, temp_file));
    int_vector<8> text;
    ASSERT_EQ(true, load_vector_from_file(text, test_file, 1));
    ::test_rank(wt, text, text.size());
}

//! Test select methods
TYPED_TEST(WtByteTest, Select)
{
    TypeParam wt;
    ASSERT_EQ(true, load_from_file(wt, temp_file));
    int_vector<8> text;
    ASSERT_EQ(true, load_vector_from_file(text, test_file, 1));
    vector<size_type> cnt(256, 0);
    ASSERT_EQ(text.size(), wt.size());
    for (size_type j=0; j<text.size(); ++j) {
        cnt[text[j]]++;
        ASSERT_EQ(j, wt.select(cnt[text[j]], text[j]))<< " j = "<<j<<" text[j]"<<text[j];
    }
}

//! Test access after swap
TYPED_TEST(WtByteTest, SwapTest)
{
    TypeParam wt1;
    ASSERT_EQ(true, load_from_file(wt1, temp_file));
    TypeParam wt2;
    wt1.swap(wt2);
    int_vector<8> text;
    ASSERT_EQ(true, load_vector_from_file(text, test_file, 1));
    ASSERT_EQ(text.size(), wt2.size());
    for (size_type j=0; j<text.size(); ++j) {
        ASSERT_EQ(wt2[j], (typename TypeParam::value_type)text[j]);
    }
}

TYPED_TEST(WtByteTest, CreatePartiallyTest)
{
    int_vector_file_buffer<8> text_buf;
    text_buf.load_from_plain(test_file);
    int_vector<8> text;
    ASSERT_EQ(true, load_vector_from_file(text, test_file, 1));
    size_type n = min(text.size(), (size_type)50);
    TypeParam wt(text_buf, n);
    ::test_rank(wt, text, n);
}

TYPED_TEST(WtByteTest, DeleteTest)
{
    std::remove(temp_file.c_str());
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " test_file temp_file" << endl;
        cout << " (1) Generates a WT out of test_file; stores it in temp_file." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
    }
    test_file = argv[1];
    temp_file  = argv[2];

    return RUN_ALL_TESTS();
}
