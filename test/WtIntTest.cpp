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

// forward declaration
template<class t_wt>
void test_interval_symbols(t_wt& wt);
// forward declaration
template<class t_wt>
void test_lex_count(t_wt& wt);
// forward declaration
template<class t_wt>
void test_lex_smaller_count(t_wt& wt);


template<class t_wt, bool lex_ordered = t_wt::lex_ordered>
struct wt_test_trait;

template<class t_wt>
struct wt_test_trait<t_wt, false> {
    static void interval_symbols_test(t_wt& wt) {
        test_interval_symbols(wt);
    }
    static void lex_count_test(t_wt&) {}
    static void lex_smaller_count_test(t_wt&) {}
};

template<class t_wt>
struct wt_test_trait<t_wt, true> {
    static void interval_symbols_test(t_wt& wt) {
        test_interval_symbols(wt);
    }
    static void lex_count_test(t_wt& wt) {
        test_lex_count(wt);
    }
    static void lex_smaller_count_test(t_wt& wt) {
        test_lex_smaller_count(wt);
    }
};

template<class t_bitvector,class t_rank, class t_select, class t_wt>
struct wt_test_trait<wt_rlmn<t_bitvector,t_rank,t_select,t_wt>,false> {
    static void interval_symbols_test(wt_rlmn<t_bitvector,t_rank,t_select,t_wt>&) {}
    static void lex_count_test(wt_rlmn<t_bitvector,t_rank,t_select,t_wt>&) {}
    static void lex_smaller_count_test(wt_rlmn<t_bitvector,t_rank,t_select,t_wt>&) {}
};


template<class T>
class WtIntTest : public ::testing::Test { };

using testing::Types;

// TODO: * add test cases for range_search_2d

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

template<class t_T>
void test_interval_symbols(t_T& wt)
{
    ASSERT_TRUE(load_from_file(wt, temp_file));

    size_type k = 0;
    std::vector<size_type> rank_c_i(wt.sigma);
    std::vector<size_type> rank_c_j(wt.sigma);
    std::vector<int_vector<>::value_type> cs(wt.sigma);

    std::mt19937_64 rng;
    for (size_type n=1; n<4; ++n) {
        std::uniform_int_distribution<uint64_t> distribution(0, n*n*n*10);
        auto dice = bind(distribution, rng);
        for (size_type i=0, j=0; i < wt.size(); i=j) {
            j = std::min(wt.size(),i+dice());

            wt.interval_symbols(i, j, k, cs, rank_c_i, rank_c_j);

            size_type symbols = (j-i);
            for (size_type m = 0; m<k; ++m) {
                ASSERT_EQ(wt.rank(i, cs[m]), rank_c_i[m]);
                ASSERT_EQ(wt.rank(j, cs[m]), rank_c_j[m]);
                ASSERT_LT((size_type)0, rank_c_j[m]-rank_c_i[m]);
                symbols -= (rank_c_j[m]-rank_c_i[m]);
                if (m>0 and t_T::lex_ordered) {
                    ASSERT_LT(cs[m-1],cs[m]);
                }
            }

            ASSERT_EQ((size_type)0, symbols);
            if (!t_T::lex_ordered) {
                sort(cs.begin(), cs.begin()+k);
                for (size_type m=1; m<k; m++) {
                    ASSERT_LT(cs[m-1], cs[m]);
                }
            }
        }
    }
}

//! Test the load method and interval_symbols method
TYPED_TEST(WtIntTest, LoadAndIntervalSymbols)
{
    TypeParam wt;
    ::wt_test_trait<TypeParam>::interval_symbols_test(wt);
}

template<class t_T>
void test_lex_count(t_T& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    std::mt19937_64 rng;
    uint64_t min = std::numeric_limits<uint64_t>::max(), max = 0;
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
            size_type smaller_c1=0, greater_c1=0, smaller_c2=0, greater_c2=0;
            int_vector<>::value_type c1=iv[i], c2=dice_symbol();
            for (; idx<j; ++idx) {
                if (iv[idx]<c1) ++smaller_c1;
                if (iv[idx]>c1) ++greater_c1;
                if (iv[idx]<c2) ++smaller_c2;
                if (iv[idx]>c2) ++greater_c2;

            }
            auto res1 = wt.lex_count(i, j, c1);
            ASSERT_EQ(wt.rank(i, c1), std::get<0>(res1));
            ASSERT_EQ(smaller_c1, std::get<1>(res1));
            ASSERT_EQ(greater_c1, std::get<2>(res1));

            auto res2 = wt.lex_count(i, j, c2);
            ASSERT_EQ(wt.rank(i, c2), std::get<0>(res2));
            ASSERT_EQ(smaller_c2, std::get<1>(res2));
            ASSERT_EQ(greater_c2, std::get<2>(res2));

            auto res3 = wt.lex_count(i, j, max+1+dice_symbol());
            ASSERT_EQ((size_type)0, std::get<0>(res3));
            ASSERT_EQ(j-i, std::get<1>(res3));
            ASSERT_EQ((size_type)0, std::get<2>(res3));
        }
    }
}

//! Test the load method and lex_count method
TYPED_TEST(WtIntTest, LoadAndLexCount)
{
    TypeParam wt;
    ::wt_test_trait<TypeParam>::lex_count_test(wt);
}

template<class t_T>
void test_lex_smaller_count(t_T& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    std::mt19937_64 rng;
    uint64_t min = std::numeric_limits<uint64_t>::max(), max = 0;
    for (size_type j=0; j < iv.size(); ++j) {
        if (min>iv[j]) min = iv[j];
        if (max<iv[j]) max = iv[j];
    }
    std::uniform_int_distribution<uint64_t> symbol_distribution(min, max);
    auto dice_symbol = bind(symbol_distribution, rng);
    int_vector<> chars(3);
    for (size_type idx=0; idx < iv.size(); ++idx) {
        chars[0] = iv[idx];
        chars[1] = dice_symbol();
        chars[2] = max+1+dice_symbol();

        for (uint64_t i = 0; i<chars.size(); ++i) {
            auto exp = wt.lex_count(0, idx, chars[i]);
            auto res = wt.lex_smaller_count(idx, chars[i]);
            ASSERT_EQ(idx-std::get<2>(exp)-std::get<1>(exp), std::get<0>(res));
            ASSERT_EQ(std::get<1>(exp), std::get<1>(res));
        }
    }
}

//! Test the load method and lex_smaller_count method
TYPED_TEST(WtIntTest, LoadAndLexSmallerCount)
{
    TypeParam wt;
    ::wt_test_trait<TypeParam>::lex_smaller_count_test(wt);
}


template<class t_wt>
void test_range_search_2d(t_wt&) {}

template<class t_bv, class t_rank, class t_sel1, class t_sel0>
void test_range_search_2d(sdsl::wt_int<t_bv, t_rank, t_sel1, t_sel0>& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);

    ASSERT_TRUE(load_from_file(wt, temp_file));

    if (wt.size() == 0)
        return;

    std::vector<uint64_t> buf(100);
    std::vector<uint64_t> unique_buf(buf.size());

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> range_distr(0, wt.size()-1);
    auto dice_range = bind(range_distr, rng);

    std::uniform_int_distribution<uint64_t> rank_distr(0, buf.size());
    auto dice_rank = bind(rank_distr, rng);

    for (size_type n=0; n<1000; ++n) {
        size_type lb = dice_range();
        size_type rb = lb+buf.size()-1;
        rb = (rb >= wt.size()) ? wt.size()-1 : rb;

        auto buf_end = std::copy(iv.begin()+lb, iv.begin()+rb+1, buf.begin());
        std::sort(buf.begin(), buf_end);
        auto unique_end = std::unique_copy(buf.begin(), buf_end,
                                           unique_buf.begin());
        size_type r1 = dice_rank() % (unique_end - unique_buf.begin());
        size_type r2 = dice_rank() % (unique_end - unique_buf.begin());
        if (r1 > r2)
            std::swap(r1, r2);
        auto vlb = unique_buf[r1];
        auto vrb = unique_buf[r2];

        size_t cnt = std::upper_bound(buf.begin(), buf_end, vrb) -
                     std::lower_bound(buf.begin(), buf_end, vlb);

        auto res = wt.range_search_2d(lb, rb, vlb, vrb);
        ASSERT_EQ(cnt, res.first);

        for (auto point : res.second) {
            // check that position is in range
            ASSERT_TRUE(point.first >= lb);
            ASSERT_TRUE(point.first <= rb);
            // check that value is in range
            ASSERT_TRUE(point.second >= vlb);
            ASSERT_TRUE(point.second <= vrb);
            // check that in the original data
            ASSERT_EQ(iv[point.first], point.second);
        }
    }
}

//! Test the load method and range_search_2d
TYPED_TEST(WtIntTest, RangeSearch2d)
{
    TypeParam wt;
    test_range_search_2d(wt);
}

template<class t_wt>
void test_quantile_freq(t_wt&) {}

template<class t_bv, class t_rank, class t_sel1, class t_sel0>
void test_quantile_freq(sdsl::wt_int<t_bv, t_rank, t_sel1, t_sel0>& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);

    ASSERT_TRUE(load_from_file(wt, temp_file));

    if (wt.size() == 0)
        return;

    std::vector<uint64_t> buf(100);

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> range_distr(0, wt.size()-1);
    auto dice_lb = bind(range_distr, rng);

    std::uniform_int_distribution<uint64_t> rank_distr(1, buf.size());
    auto dice_range = bind(rank_distr, rng);

    for (size_type n=0; n<1000; ++n) {
        size_type lb = dice_lb();
        size_type rb = lb+dice_range()-1;
        rb = (rb >= wt.size()) ? wt.size()-1 : rb;

        auto buf_end = std::copy(iv.begin()+lb, iv.begin()+rb+1, buf.begin());
        std::sort(buf.begin(), buf_end);

        for (auto it = buf.begin(); it!=buf_end; ++it) {
            auto val = *it;
            size_type freq = std::upper_bound(buf.begin(), buf_end, val) -
                             std::lower_bound(buf.begin(), buf_end, val);
            size_type q    = it - buf.begin();
            auto res = quantile_freq(wt, lb, rb, q);
            ASSERT_EQ(val, res.first);
            ASSERT_EQ(freq, res.second);
        }
    }
}

//! Test the load method and quantile_freq
TYPED_TEST(WtIntTest, QuantileFreq)
{
    TypeParam wt;
    test_quantile_freq(wt);
}

template<class t_wt>
void test_intersect(t_wt&) {}

template<class t_bv, class t_rank, class t_sel1, class t_sel0>
void test_intersect(sdsl::wt_int<t_bv, t_rank, t_sel1, t_sel0>& wt)
{
    typedef sdsl::wt_int<t_bv, t_rank, t_sel1, t_sel0> t_wt;
    using t_pvs = std::pair<typename t_wt::value_type,
          typename t_wt::size_type>;
    int_vector<> iv;
    load_from_file(iv, test_file);

    ASSERT_TRUE(load_from_file(wt, temp_file));

    if (wt.size() == 0)
        return;

    std::vector<std::vector<uint64_t>> buf(2, std::vector<uint64_t>(300));
    std::vector<uint64_t> res_buf(buf[0].size()*2);
    std::vector<std::vector<uint64_t>::iterator> buf_end(2);

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> range_distr(0, wt.size()-1);
    auto dice_lb = bind(range_distr, rng);

    std::uniform_int_distribution<uint64_t> rank_distr(1, buf[0].size()-1);
    auto dice_range = bind(rank_distr, rng);

    for (size_type n=0; n<1000; ++n) {
        range_vec_type ranges;

        for (size_t i=0; i<2; ++i) {
            size_type lb = dice_lb();
            size_type rb = lb+dice_range()-1;
            rb = (rb >= wt.size()) ? wt.size()-1 : rb;
            ranges.emplace_back(lb,rb);
            buf_end[i] = std::copy(iv.begin()+lb,
                                   iv.begin()+rb+1, buf[i].begin());
            sort(buf[i].begin(), buf_end[i]);
        }

        auto res_end =
            std::set_intersection(buf[0].begin(), buf_end[0],
                                  buf[1].begin(), buf_end[1],
                                  res_buf.begin());
        res_end = std::unique(res_buf.begin(), res_end);

        auto itsct = intersect(wt, ranges);
        size_type res_size = res_end-res_buf.begin();
        ASSERT_EQ(res_size, itsct.size());
        if (!t_wt::lex_ordered) {
            sort(itsct.begin(), itsct.end(), [](const t_pvs& x, const t_pvs& y) {
                return x.first < y.first;
            });
        }
        for (size_t i=0; i< itsct.size(); ++i) {
            ASSERT_EQ(res_buf[i], itsct[i].first);
        }
    }
}

//! Test the load method and intersect
TYPED_TEST(WtIntTest, Intersect)
{
    TypeParam wt;
    test_intersect(wt);
}

TYPED_TEST(WtIntTest, DeleteTest)
{
    sdsl::remove(temp_file);
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
