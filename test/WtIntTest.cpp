#include "sdsl/wavelet_trees.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <map>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;
typedef std::map<int_vector<>::value_type,size_type> tMII;

string test_file;
string temp_file;
bool in_memory;

template<class T>
class WtIntTest : public ::testing::Test { };

using testing::Types;

// TODO: * add test cases for range_search_2d
//       * add test for lex_count,...

typedef Types<
wt_blcd<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_huff<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_huff<rrr_vector<63>, rrr_vector<63>::rank_1_type, rrr_vector<63>::select_1_type, rrr_vector<63>::select_0_type, int_tree<>>
        ,wt_hutu<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_int<rrr_vector<15> >
        ,wt_int<>
        ,wt_int<rrr_vector<63> >
        > Implementations;

TYPED_TEST_CASE(WtIntTest, Implementations);

//! Test the parametrized constructor
TYPED_TEST(WtIntTest, Constructor)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    double iv_size = size_in_mega_bytes(iv);
    std::cout << "tc = " << test_file << std::endl;
    {
        TypeParam wt;
        sdsl::construct(wt, test_file);
        std::cout << "compression = " << size_in_mega_bytes(wt)/iv_size << std::endl;
        ASSERT_EQ(iv.size(), wt.size());
        for (size_type j=0; j < iv.size(); ++j) {
            ASSERT_EQ(iv[j], wt[j])<<j;
        }
        ASSERT_TRUE(store_to_file(wt, temp_file));
    }
    {
        int_vector_buffer<> iv_buf(test_file);
        TypeParam wt(iv_buf, 0);
        ASSERT_EQ((size_type)0,  wt.size());
    }
    {
        int_vector_buffer<> iv_buf(test_file);
        size_type len = (iv.size() >= 6) ? 6 : iv.size();
        TypeParam wt(iv_buf, len);
        ASSERT_EQ(len, wt.size());
        for (size_type j=0; j < len; ++j) {
            ASSERT_EQ(iv[j], wt[j])<<j;
        }
    }
}

//! Test loading and accessing the wavelet tree
TYPED_TEST(WtIntTest, LoadAndAccess)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    for (size_type j=0; j < iv.size(); ++j) {
        ASSERT_EQ(iv[j], wt[j])<<j;
    }
}

//! Test the load method and rank method
TYPED_TEST(WtIntTest, LoadAndRank)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    tMII check_rank;
    for (size_type j=0; j < iv.size(); ++j) {
        ASSERT_EQ(wt.rank(j, iv[j]), check_rank[iv[j]]);
        check_rank[iv[j]]++;
    }
    for (tMII::const_iterator it=check_rank.begin(); it!=check_rank.end(); ++it) {
        ASSERT_EQ(wt.rank(wt.size(), it->first), it->second);
    }
}

//! Test the load method and select method
TYPED_TEST(WtIntTest, LoadAndSelect)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    tMII count;
    for (size_type j=0; j < iv.size(); ++j) {
        count[iv[j]]++;
        ASSERT_EQ(j, wt.select(count[iv[j]], iv[j]))
                << "iv[j]=" << iv[j] << " j="<<j;
    }
}

//! Test the load method and inverse_select method
TYPED_TEST(WtIntTest, LoadAndInverseSelect)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    tMII check_rank;
    for (size_type j=0; j < iv.size(); ++j) {
        typename TypeParam::value_type x;
        ASSERT_EQ(check_rank[iv[j]], wt.inverse_select(j, x));
        ASSERT_EQ(iv[j], x);
        check_rank[iv[j]]++;
    }
}

TYPED_TEST(WtIntTest, DeleteTest)
{
    sdsl::remove(temp_file);
}

template<class T>
class WtIntTopK : public ::testing::Test { };
typedef Types< wt_int<rrr_vector<15>> , wt_int<> , wt_int<rrr_vector<63>>  > Implementations_ordered;

TYPED_TEST_CASE(WtIntTopK, Implementations_ordered);

TYPED_TEST(WtIntTopK, quantile_freq)
{
    int_vector<> a = {7,2,2,4,4,4,15,14,17,3,8,10,11,1,18,5,6,4,16,9,12,13};

    temp_file = ram_file_name("test.file");
    store_to_file(a, temp_file);

    TypeParam wt;
    sdsl::construct(wt, temp_file);

    auto res = wt.quantile_freq(0,8,0);
    ASSERT_EQ((size_type)2,res.first);
    ASSERT_EQ((size_type)2,res.second);
    res = wt.quantile_freq(0,8,1);
    ASSERT_EQ((size_type)2, res.first);
    ASSERT_EQ((size_type)2, res.second);
    res = wt.quantile_freq(0,8,2);
    ASSERT_EQ((size_type)4, res.first);
    ASSERT_EQ((size_type)3, res.second);
    res = wt.quantile_freq(0,8,3);
    ASSERT_EQ((size_type)4, res.first);
    ASSERT_EQ((size_type)3, res.second);
    res = wt.quantile_freq(0,8,4);
    ASSERT_EQ((size_type)4, res.first);
    ASSERT_EQ((size_type)3, res.second);
    res = wt.quantile_freq(0,8,5);
    ASSERT_EQ((size_type)7,res.first);
    ASSERT_EQ((size_type)1,res.second);
    res = wt.quantile_freq(0,8,6);
    ASSERT_EQ((size_type)14,res.first);
    ASSERT_EQ((size_type)1,res.second);
    res = wt.quantile_freq(0,8,7);
    ASSERT_EQ((size_type)15,res.first);
    ASSERT_EQ((size_type)1,res.second);
    res = wt.quantile_freq(0,8,8);
    ASSERT_EQ((size_type)17,res.first);
    ASSERT_EQ((size_type)1,res.second);

}


TYPED_TEST(WtIntTopK, topk_greedy)
{
    int_vector<> a = {2,2,2,4,4,4,3,4,17,3,8,10,11,1,18,5,6,4,16,9,12,13,7,14,15,0};

    auto f = ram_file_name("test.file");
    store_to_file(a, f);

    TypeParam wt;
    sdsl::construct(wt, f);

    auto results = wt.topk_greedy(0,9,3);

    ASSERT_EQ((size_type)3,results.size());
    ASSERT_EQ((size_type)4,results[0].first);
    ASSERT_EQ((size_type)4,results[0].second);
    ASSERT_EQ((size_type)2,results[1].first);
    ASSERT_EQ((size_type)3,results[1].second);
    ASSERT_EQ((size_type)3,results[2].first);
    ASSERT_EQ((size_type)2,results[2].second);

}

TYPED_TEST(WtIntTopK, topk_qprobing)
{
    int_vector<> a = {2,2,2,4,4,4,3,4,17,3,8,10,11,1,18,5,6,4,16,9,12,13,7,14,15,0};

    auto f = ram_file_name("test.file");
    store_to_file(a, f);

    TypeParam wt;
    sdsl::construct(wt, f);

    auto results = wt.topk_qprobing(0,9,3);

    ASSERT_EQ((size_type)3,results.size());
    ASSERT_EQ((size_type)4,results[0].first);
    ASSERT_EQ((size_type)4,results[0].second);
    ASSERT_EQ((size_type)2,results[1].first);
    ASSERT_EQ((size_type)3,results[1].second);
    ASSERT_EQ((size_type)3,results[2].first);
    ASSERT_EQ((size_type)2,results[2].second);

    results = wt.topk_qprobing(0,9,5);

    ASSERT_EQ((size_type)4,results.size());
    ASSERT_EQ((size_type)4,results[0].first);
    ASSERT_EQ((size_type)4,results[0].second);
    ASSERT_EQ((size_type)2,results[1].first);
    ASSERT_EQ((size_type)3,results[1].second);
    ASSERT_EQ((size_type)3,results[2].first);
    ASSERT_EQ((size_type)2,results[2].second);
    ASSERT_EQ((size_type)17,results[3].first);
    ASSERT_EQ((size_type)1,results[3].second);
}

TYPED_TEST(WtIntTopK, intersection)
{
    int_vector<> a = {2,2,2,4,5,4,0,0,4,5,7,7,6,6,3,5,1,6};

    auto f = ram_file_name("test.file");
    store_to_file(a, f);

    TypeParam wt;
    sdsl::construct(wt, f);

    std::vector< std::pair<size_type,size_type> > ranges;
    ranges.push_back(std::pair<size_type,size_type>(2,4));
    ranges.push_back(std::pair<size_type,size_type>(6,9));

    auto results = wt.intersect(ranges);

    ASSERT_EQ((size_type)2, results.size());
    pair<size_type,size_type> p1 = results[0];
    ASSERT_EQ((size_type)5, p1.first);
    pair<size_type,size_type> p2 = results[1];
    ASSERT_EQ((size_type)4, p2.first);
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " test_file temp_file [in-memory]" << endl;
        cout << " (1) Generates a WT out of test_file; stores it in temp_file." << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
    }
    test_file = argv[1];
    temp_file = argv[2];
    in_memory = argc > 3;
    if (in_memory) {
        int_vector<> data;
        load_from_file(data, test_file);
        test_file = ram_file_name(test_file);
        store_to_file(data, test_file);
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}
