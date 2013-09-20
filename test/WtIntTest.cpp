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
//       * add test for lex_smaller_count,...

typedef Types<
wt_blcd<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_huff<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_huff<rrr_vector<63>, rrr_vector<63>::rank_1_type, rrr_vector<63>::select_1_type, rrr_vector<63>::select_0_type, int_tree<>>
        ,wt_hutu<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_int<rrr_vector<15>>
        ,wt_int<>
        ,wt_int<rrr_vector<63>>
        ,wt_rlmn<bit_vector, rank_support_v5<>, select_support_mcl<1>, wt_int<>>
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
        std::set<uint64_t> sigma_set;
        for (size_type j=0; j < iv.size(); ++j) {
            ASSERT_EQ(iv[j], wt[j])<<j;
            sigma_set.insert(iv[j]);
        }
        ASSERT_EQ(sigma_set.size(), wt.sigma);
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
        auto rc = wt.inverse_select(j);
        ASSERT_EQ(check_rank[iv[j]], rc.first);
        ASSERT_EQ(iv[j], rc.second);
        check_rank[iv[j]]++;
    }
}


TYPED_TEST(WtIntTest, DeleteTest)
{
    sdsl::remove(temp_file);
}

template<class T>
class WtIntLexOrdered : public ::testing::Test { };
typedef Types<
wt_blcd<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_hutu<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_int<>
        ,wt_int<rrr_vector<15>>
        ,wt_int<rrr_vector<63>>
        > Implementations_lex_ordered;

TYPED_TEST_CASE(WtIntLexOrdered, Implementations_lex_ordered);

//! Test the parametrized constructor
TYPED_TEST(WtIntLexOrdered, Constructor)
{
    TypeParam wt;
    sdsl::construct(wt, test_file);
    ASSERT_TRUE(store_to_file(wt, temp_file));
}

//! Test the load method and lex_count method
TYPED_TEST(WtIntLexOrdered, LoadAndLexCount)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    TypeParam wt;
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    std::mt19937_64 rng;
    uint64_t min = UINT64_MAX, max = 0;
    for (size_type j=0; j < iv.size(); ++j) {
        if (min>iv[j]) min = iv[j];
        if (max<iv[j]) max = iv[j];
    }
    std::uniform_int_distribution<uint64_t> symbol_distribution(min, max);
    auto dice_symbol = bind(symbol_distribution, rng);
    for (size_type k=1; k<4; ++k) {
        std::uniform_int_distribution<uint64_t> distribution(0, k*k*k*10);
        auto dice = bind(distribution, rng);
        for (size_type idx=0; idx < iv.size();) {
            size_type i = idx, j = std::min(wt.size(),i+dice());
            size_type smaller_c1=0,greater_c1=0,smaller_c2=0,greater_c2=0;
            int_vector<>::value_type c1=iv[i],c2=dice_symbol();
            for (; idx<j; ++idx) {
                if (iv[idx]<c1) ++smaller_c1;
                if (iv[idx]>c1) ++greater_c1;
                if (iv[idx]<c2) ++smaller_c2;
                if (iv[idx]>c2) ++greater_c2;

            }
            auto res1 = wt.lex_count(i,j,c1);
            ASSERT_EQ(wt.rank(i,c1),std::get<0>(res1));
            ASSERT_EQ(smaller_c1,std::get<1>(res1));
            ASSERT_EQ(greater_c1,std::get<2>(res1));

            auto res2 = wt.lex_count(i,j,c2);
            ASSERT_EQ(wt.rank(i,c2),std::get<0>(res2));
            ASSERT_EQ(smaller_c2,std::get<1>(res2));
            ASSERT_EQ(greater_c2,std::get<2>(res2));

            auto res3 = wt.lex_count(i,j,max+1);
            ASSERT_EQ(0,std::get<0>(res3));
            ASSERT_EQ(j-i,std::get<1>(res3));
            ASSERT_EQ(0,std::get<2>(res3));
        }
    }
}

TYPED_TEST(WtIntLexOrdered, DeleteTest)
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
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " test_file temp_file [in-memory]" << endl;
        cout << " (1) Generates a WT out of test_file; stores it in temp_file." << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
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
