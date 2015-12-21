#include "sdsl/wavelet_trees.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <map>
#include <queue>
#include <algorithm>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;
typedef map<int_vector<>::value_type,size_type> tMII;

string test_file;
string temp_file;
string temp_dir;
bool in_memory;

template<class T>
class wt_int_test : public ::testing::Test { };

using testing::Types;

typedef Types<
wt_blcd<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_huff<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_huff<rrr_vector<63>, rrr_vector<63>::rank_1_type, rrr_vector<63>::select_1_type, rrr_vector<63>::select_0_type, int_tree<>>
        ,wt_hutu<bit_vector, rank_support_v<>, select_support_mcl<1>, select_support_mcl<0>, int_tree<>>
        ,wt_gmr<>
        ,wt_ap<>
        ,wm_int<>
        ,wt_int<>
        ,wt_int<rrr_vector<15>>
        ,wt_int<rrr_vector<63>>
        ,wt_rlmn<bit_vector, rank_support_v5<>, select_support_mcl<1>, wt_int<>>
        > Implementations;

TYPED_TEST_CASE(wt_int_test, Implementations);

//! Test the parametrized constructor
TYPED_TEST(wt_int_test, constructor)
{
    static_assert(sdsl::util::is_regular<TypeParam>::value, "Type is not regular");
    int_vector<> iv;
    load_from_file(iv, test_file);
    double iv_size = size_in_mega_bytes(iv);
    cout << "tc = " << test_file << endl;
    {
        TypeParam wt;
        sdsl::construct(wt, test_file);
        cout << "compression = " << size_in_mega_bytes(wt)/iv_size << endl;
        ASSERT_EQ(iv.size(), wt.size());
        set<uint64_t> sigma_set;
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
TYPED_TEST(wt_int_test, load_and_access)
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
TYPED_TEST(wt_int_test, load_and_rank)
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
    for (auto it=check_rank.begin(); it!=check_rank.end(); ++it) {
        ASSERT_EQ(wt.rank(wt.size(), it->first), it->second);
    }
}

//! Test the load method and rank method
TYPED_TEST(wt_int_test, load_and_move_and_rank)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    TypeParam wt_load;
    ASSERT_TRUE(load_from_file(wt_load, temp_file));
    TypeParam wt = move(wt_load);
    ASSERT_EQ(iv.size(), wt.size());
    tMII check_rank;
    for (size_type j=0; j < iv.size(); ++j) {
        ASSERT_EQ(wt.rank(j, iv[j]), check_rank[iv[j]]);
        check_rank[iv[j]]++;
    }
    for (auto it=check_rank.begin(); it!=check_rank.end(); ++it) {
        ASSERT_EQ(wt.rank(wt.size(), it->first), it->second);
    }
}

//! Test the load method and select method
TYPED_TEST(wt_int_test, load_and_select)
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

//! Test the load method and select method
TYPED_TEST(wt_int_test, load_and_move_and_select)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    TypeParam wt_load;
    ASSERT_TRUE(load_from_file(wt_load, temp_file));
    TypeParam wt = move(wt_load);
    ASSERT_EQ(iv.size(), wt.size());
    tMII count;
    for (size_type j=0; j < iv.size(); ++j) {
        count[iv[j]]++;
        ASSERT_EQ(j, wt.select(count[iv[j]], iv[j]))
                << "iv[j]=" << iv[j] << " j="<<j;
    }
}

//! Test the load method and inverse_select method
TYPED_TEST(wt_int_test, load_and_inverse_select)
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

template<class t_wt>
void
test_interval_symbols(typename enable_if<!(has_node_type<t_wt>::value),
                      t_wt>::type&)
{
    // interval_symbols not implemented
}

template<class t_wt>
void
test_interval_symbols(typename enable_if<has_node_type<t_wt>::value,
                      t_wt>::type& wt)
{
    ASSERT_TRUE(load_from_file(wt, temp_file));

    size_type k = 0;
    vector<size_type> rank_c_i(wt.sigma);
    vector<size_type> rank_c_j(wt.sigma);
    vector<int_vector<>::value_type> cs(wt.sigma);

    mt19937_64 rng;
    for (size_type n=1; n<4; ++n) {
        uniform_int_distribution<uint64_t> distribution(0, n*n*n*10);
        auto dice = bind(distribution, rng);
        for (size_type i=0, j=0; i < wt.size(); i=j) {
            j = min(wt.size(),i+dice());

            interval_symbols(wt, i, j, k, cs, rank_c_i, rank_c_j);

            size_type symbols = (j-i);
            for (size_type m = 0; m<k; ++m) {
                ASSERT_EQ(wt.rank(i, cs[m]), rank_c_i[m]);
                ASSERT_EQ(wt.rank(j, cs[m]), rank_c_j[m]);
                ASSERT_LT((size_type)0, rank_c_j[m]-rank_c_i[m]);
                symbols -= (rank_c_j[m]-rank_c_i[m]);
                if (m>0 and t_wt::lex_ordered) {
                    ASSERT_LT(cs[m-1],cs[m]);
                }
            }

            ASSERT_EQ((size_type)0, symbols);
            if (!t_wt::lex_ordered) {
                sort(cs.begin(), cs.begin()+k);
                for (size_type m=1; m<k; m++) {
                    ASSERT_LT(cs[m-1], cs[m]);
                }
            }
        }
    }
}

//! Test the load method and interval_symbols method
TYPED_TEST(wt_int_test, load_and_interval_symbols)
{
    TypeParam wt;
    test_interval_symbols<TypeParam>(wt);
}

template<class t_wt>
void
test_lex_count(typename enable_if<!(t_wt::lex_ordered), t_wt>::type&)
{
    // lex_count not implemented
}

template<class t_wt>
void
test_lex_count(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    mt19937_64 rng;
    uint64_t min = numeric_limits<uint64_t>::max(), max = 0;
    for (size_type j=0; j < iv.size(); ++j) {
        if (min>iv[j]) min = iv[j];
        if (max<iv[j]) max = iv[j];
    }
    uniform_int_distribution<uint64_t> symbol_distribution(min, max);
    auto dice_symbol = bind(symbol_distribution, rng);
    for (size_type k=1; k<4; ++k) {
        uniform_int_distribution<uint64_t> distribution(0, k*k*k*10);
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
            ASSERT_EQ(wt.rank(i, c1), get<0>(res1));
            ASSERT_EQ(smaller_c1, get<1>(res1));
            ASSERT_EQ(greater_c1, get<2>(res1));

            auto res2 = wt.lex_count(i, j, c2);
            ASSERT_EQ(wt.rank(i, c2), get<0>(res2));
            ASSERT_EQ(smaller_c2, get<1>(res2));
            ASSERT_EQ(greater_c2, get<2>(res2));

            auto res3 = wt.lex_count(i, j, max+1+dice_symbol());
            ASSERT_EQ((size_type)0, get<0>(res3));
            ASSERT_EQ(j-i, get<1>(res3));
            ASSERT_EQ((size_type)0, get<2>(res3));
        }
    }
}

//! Test the load method and lex_count method
TYPED_TEST(wt_int_test, load_and_lex_count)
{
    TypeParam wt;
    test_lex_count<TypeParam>(wt);
}

template<class t_wt>
void
test_lex_smaller_count(typename enable_if<!(t_wt::lex_ordered), t_wt>::type&)
{
    // lex_smaller_count not implemented
}

template<class t_wt>
void
test_lex_smaller_count(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    mt19937_64 rng;
    uint64_t min = numeric_limits<uint64_t>::max(), max = 0;
    for (size_type j=0; j < iv.size(); ++j) {
        if (min>iv[j]) min = iv[j];
        if (max<iv[j]) max = iv[j];
    }
    uniform_int_distribution<uint64_t> symbol_distribution(min, max);
    auto dice_symbol = bind(symbol_distribution, rng);
    int_vector<> chars(3);
    for (size_type idx=0; idx < iv.size(); ++idx) {
        chars[0] = iv[idx];
        chars[1] = dice_symbol();
        chars[2] = max+1+dice_symbol();

        for (uint64_t i = 0; i<chars.size(); ++i) {
            auto exp = wt.lex_count(0, idx, chars[i]);
            auto res = wt.lex_smaller_count(idx, chars[i]);
            ASSERT_EQ(idx-get<2>(exp)-get<1>(exp), get<0>(res));
            ASSERT_EQ(get<1>(exp), get<1>(res));
        }
    }
}

//! Test the load method and lex_smaller_count method
TYPED_TEST(wt_int_test, load_and_lex_smaller_count)
{
    TypeParam wt;
    test_lex_smaller_count<TypeParam>(wt);
}


template<class t_wt>
void
test_range_search_2d(typename enable_if<!(has_range_search_2d<t_wt>::value),
                     t_wt>::type&) {}

template<class t_wt>
void
test_range_search_2d(typename enable_if<has_range_search_2d<t_wt>::value,
                     t_wt>::type& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);

    ASSERT_TRUE(load_from_file(wt, temp_file));

    if (wt.size() == 0)
        return;

    vector<uint64_t> buf(100);
    vector<uint64_t> unique_buf(buf.size());

    mt19937_64 rng;
    uniform_int_distribution<uint64_t> range_distr(0, wt.size()-1);
    auto dice_range = bind(range_distr, rng);

    uniform_int_distribution<uint64_t> rank_distr(0, buf.size());
    auto dice_rank = bind(rank_distr, rng);

    for (size_type n=0; n<1000; ++n) {
        size_type lb = dice_range();
        size_type rb = lb+buf.size()-1;
        rb = (rb >= wt.size()) ? wt.size()-1 : rb;

        auto buf_end = copy(iv.begin()+lb, iv.begin()+rb+1, buf.begin());
        sort(buf.begin(), buf_end);
        auto unique_end = unique_copy(buf.begin(), buf_end,
                                      unique_buf.begin());
        size_type r1 = dice_rank() % (unique_end - unique_buf.begin());
        size_type r2 = dice_rank() % (unique_end - unique_buf.begin());
        if (r1 > r2)
            swap(r1, r2);
        auto vlb = unique_buf[r1];
        auto vrb = unique_buf[r2];

        size_t cnt = upper_bound(buf.begin(), buf_end, vrb) -
                     lower_bound(buf.begin(), buf_end, vlb);

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
TYPED_TEST(wt_int_test, range_search_2d)
{
    TypeParam wt;
    test_range_search_2d<TypeParam>(wt);
}

template<class t_wt>
void
test_quantile_freq(typename enable_if<!t_wt::lex_ordered, t_wt>::type&) {}

template<class t_wt>
void
test_quantile_freq(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);

    ASSERT_TRUE(load_from_file(wt, temp_file));

    if (wt.size() == 0)
        return;

    vector<uint64_t> buf(100);

    mt19937_64 rng;
    uniform_int_distribution<uint64_t> range_distr(0, wt.size()-1);
    auto dice_lb = bind(range_distr, rng);

    uniform_int_distribution<uint64_t> rank_distr(1, buf.size());
    auto dice_range = bind(rank_distr, rng);

    for (size_type n=0; n<1000; ++n) {
        size_type lb = dice_lb();
        size_type rb = lb+dice_range()-1;
        rb = (rb >= wt.size()) ? wt.size()-1 : rb;

        auto buf_end = copy(iv.begin()+lb, iv.begin()+rb+1, buf.begin());
        sort(buf.begin(), buf_end);

        for (auto it = buf.begin(); it!=buf_end; ++it) {
            auto val = *it;
            size_type freq = upper_bound(buf.begin(), buf_end, val) -
                             lower_bound(buf.begin(), buf_end, val);
            size_type q    = it - buf.begin();
            auto res = quantile_freq(wt, lb, rb, q);
            ASSERT_EQ(val, res.first);
            ASSERT_EQ(freq, res.second);
        }
    }
}

//! Test the load method and quantile_freq
TYPED_TEST(wt_int_test, quantile_freq)
{
    TypeParam wt;
    test_quantile_freq<TypeParam>(wt);
}


template<class t_wt>
void
test_intersect(typename enable_if<!(has_node_type<t_wt>::value),t_wt>::type&)
{
    // intersect not implemented
}

template<class t_wt>
void
test_intersect(typename enable_if<has_node_type<t_wt>::value,t_wt>::type& wt)
{
    using t_pvs = pair<typename t_wt::value_type,
          typename t_wt::size_type>;
    int_vector<> iv;
    load_from_file(iv, test_file);

    ASSERT_TRUE(load_from_file(wt, temp_file));

    if (wt.size() == 0)
        return;

    vector<vector<uint64_t>> buf(2, vector<uint64_t>(300));
    vector<uint64_t> res_buf(buf[0].size()*2);
    vector<vector<uint64_t>::iterator> buf_end(2);

    mt19937_64 rng;
    uniform_int_distribution<uint64_t> range_distr(0, wt.size()-1);
    auto dice_lb = bind(range_distr, rng);

    uniform_int_distribution<uint64_t> rank_distr(1, buf[0].size()-1);
    auto dice_range = bind(rank_distr, rng);

    for (size_type n=0; n<1000; ++n) {
        range_vec_type ranges;

        for (size_t i=0; i<2; ++i) {
            size_type lb = dice_lb();
            size_type rb = lb+dice_range()-1;
            rb = (rb >= wt.size()) ? wt.size()-1 : rb;
            ranges.push_back({lb,rb});
            buf_end[i] = copy(iv.begin()+lb,
                              iv.begin()+rb+1, buf[i].begin());
            sort(buf[i].begin(), buf_end[i]);
        }

        auto res_end =
            set_intersection(buf[0].begin(), buf_end[0],
                             buf[1].begin(), buf_end[1],
                             res_buf.begin());
        res_end = unique(res_buf.begin(), res_end);

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
TYPED_TEST(wt_int_test, intersect)
{
    TypeParam wt;
    test_intersect<TypeParam>(wt);
}

template<class t_wt>
void
test_nodes(typename enable_if<!(has_node_type<t_wt>::value),t_wt>::type&)
{
    // not implemented
}

template<class t_wt>
void
test_nodes(typename enable_if<has_node_type<t_wt>::value,t_wt>::type& wt)
{
    using node_type = typename t_wt::node_type;
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    if (wt.size() < 1)
        return;
    std::vector<node_type> v;
    {
        std::queue<node_type> q;
        q.push(wt.root());
        while (!q.empty() and v.size() < 100000) {
            node_type x = q.front();
            q.pop();
            if (!wt.is_leaf(x)) {
                auto children = wt.expand(x);
                if (!wt.empty(std::get<0>(children)))
                    q.emplace(std::get<0>(children));
                if (!wt.empty(std::get<1>(children)))
                    q.emplace(std::get<1>(children));
            }
            v.push_back(x);
        }
    }
    std::sort(v.begin(), v.end());
    for (size_t i=1; i<v.size(); ++i) {
        ASSERT_TRUE(v[i-1] < v[i]);
    }
}

//! Test the load method and nodes
TYPED_TEST(wt_int_test, nodes)
{
    TypeParam wt;
    test_nodes<TypeParam>(wt);
}

template<class t_wt>
void
test_symbol_gte(typename enable_if<!(t_wt::lex_ordered), t_wt>::type&)
{
    // symbol_gte not implemented
}

template<class t_wt>
void
test_symbol_gte(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    mt19937_64 rng;
    uint64_t min = numeric_limits<uint64_t>::max(), max = 0;
    std::set<uint64_t> syms;
    for (size_type j=0; j < iv.size(); ++j) {
        if (min>iv[j]) min = iv[j];
        if (max<iv[j]) max = iv[j];
        syms.insert(iv[j]);
    }

    if (iv.size() == 0) {
        return;
    }

    // check symbols that are in there also are reported as "equal"
    auto itr = syms.begin();
    auto end = syms.end();
    while (itr != end) {
        auto value = *itr;
        auto ret = symbol_gte(wt,value);
        ASSERT_TRUE(ret.first);
        ASSERT_EQ(value,ret.second);
        ++itr;
    }

    // check symbols symbols that are smaller than than min
    for (size_t i=0; i<min; i++) {
        auto ret = symbol_gte(wt,i);
        ASSERT_TRUE(ret.first);
        ASSERT_EQ(ret.second,min);
    }

    // check symbols that are largest than max
    for (size_t i=max+100; i>max; i--) {
        auto ret = symbol_gte(wt,i);
        ASSERT_FALSE(ret.first);
    }

    // check values in between that do not exist
    for (size_t i=min; i<max; i++) {
        auto itr = syms.find(i);
        if (itr == syms.end()) {
            size_t j=i+1;
            auto next = syms.find(j);
            while (next == syms.end()) {
                next = syms.find(j+1);
                j++;
            }
            if (next != syms.end()) {
                auto next_val = *next;
                auto ret = symbol_gte(wt,i);
                ASSERT_TRUE(ret.first);
                ASSERT_EQ(ret.second,next_val);
            }
        }
    }
}

//! Test the load method and symbol_gte
TYPED_TEST(wt_int_test, symbol_gte)
{
    TypeParam wt;
    test_symbol_gte<TypeParam>(wt);
}

template<class t_wt>
void
test_symbol_lte(typename enable_if<!(t_wt::lex_ordered), t_wt>::type&)
{
    // symbol_lte not implemented
}



template<class t_wt>
void
test_symbol_lte(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt)
{
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    mt19937_64 rng;
    uint64_t min = numeric_limits<uint64_t>::max(), max = 0;
    std::set<uint64_t> syms;
    for (size_type j=0; j < iv.size(); ++j) {
        if (min>iv[j]) min = iv[j];
        if (max<iv[j]) max = iv[j];
        syms.insert(iv[j]);
    }

    if (iv.size() == 0) {
        return;
    }

    // check symbols that are in there also are reported as "equal"
    auto itr = syms.begin();
    auto end = syms.end();
    while (itr != end) {
        auto value = *itr;
        auto ret = symbol_lte(wt,value);
        ASSERT_TRUE(ret.first);
        ASSERT_EQ(value,ret.second);
        ++itr;
    }

    // check symbols symbols that are smaller than than min
    for (size_t i=0; i<min; i++) {
        auto ret = symbol_lte(wt,i);
        ASSERT_FALSE(ret.first);
        //ASSERT_EQ(ret.second,min);
    }

    // check symbols that are larget than max
    for (size_t i=max+100; i>max; i--) {
        auto ret = symbol_lte(wt,i);
        ASSERT_TRUE(ret.first);
        ASSERT_EQ(ret.second,max);
    }

    // check values in between that do not exist
    for (size_t i=min+1; i<max; i++) {
        auto itr = syms.find(i);
        if (itr == syms.end()) {
            size_t j=i-1;
            auto prev = syms.find(j);
            while (prev == syms.end()) {
                prev = syms.find(j-1);
                j--;
            }
            if (prev != syms.end()) {
                auto prev_val = *prev;
                auto ret = symbol_lte(wt,i);
                ASSERT_TRUE(ret.first);
                ASSERT_EQ(ret.second,prev_val);
            }
        }
    }
}


//! Test the load method and symbol_lte
TYPED_TEST(wt_int_test, symbol_lte)
{
    TypeParam wt;
    test_symbol_lte<TypeParam>(wt);
}


template<class t_wt>
void
test_range_unique_values(typename enable_if<!(t_wt::lex_ordered), t_wt>::type&)
{
    // test_range_unique_values not implemented
}

template<class t_wt>
void
test_range_unique_values(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt)
{
    using value_type = typename t_wt::value_type;
    int_vector<> iv;
    load_from_file(iv, test_file);
    ASSERT_TRUE(load_from_file(wt, temp_file));
    ASSERT_EQ(iv.size(), wt.size());
    value_type min = numeric_limits<value_type>::max(), max = 0;
    std::set<value_type> syms;
    for (size_type j=0; j < iv.size(); ++j) {
        if (min>iv[j]) min = iv[j];
        if (max<iv[j]) max = iv[j];
        syms.insert(iv[j]);
    }

    if (iv.size() == 0) {
        return;
    }

    // try 128 random queries
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> x_dist(0, wt.size()-1);
    std::uniform_int_distribution<uint64_t> y_dist(0, max);
    auto xdice = bind(x_dist, rng);
    auto ydice = bind(y_dist, rng);
    for (size_t i=0; i<128; i++) {
        size_t x_i = xdice();
        size_t x_j = xdice();
        if (x_i>x_j) std::swap(x_i,x_j);
        size_t y_i = ydice();
        size_t y_j = ydice();
        if (y_i>y_j) std::swap(y_i,y_j);
        auto uniq_values = restricted_unique_range_values(wt,x_i,x_j,y_i,y_j);

        /* verify */
        std::set<value_type> syms;
        for (size_t j=x_i; j<=x_j; j++) {
            if (iv[j] >= y_i && iv[j] <= y_j) syms.insert(iv[j]);
        }
        auto itr = syms.begin();
        auto end = syms.end();
        size_t r = 0;
        while (itr != end) {
            auto value = *itr;
            ASSERT_EQ(value,uniq_values[r]);
            r++;
            itr++;
        }
    }

    // check invalid queries don't do the wrong thing
    //                0 1 2 3 4 5  6  7  8  9  10
    int_vector<> S = {5,6,7,8,9,5,6,7,13,14,15};
    t_wt wt_test;
    construct_im(wt_test,S);
    auto empty_uniq_values = restricted_unique_range_values(wt,8,3,7,9);
    ASSERT_TRUE(empty_uniq_values.size() == 0);
    empty_uniq_values = restricted_unique_range_values(wt,3,8,9,7);
    ASSERT_TRUE(empty_uniq_values.size() == 0);
    empty_uniq_values = restricted_unique_range_values(wt,3,8,3,4);
    ASSERT_TRUE(empty_uniq_values.size() == 0);
    empty_uniq_values = restricted_unique_range_values(wt,3,8,10,11);
    ASSERT_TRUE(empty_uniq_values.size() == 0);
}

//! Test the load method and restricted_unique_range_values
TYPED_TEST(wt_int_test, restricted_unique_range_values)
{
    TypeParam wt;
    test_range_unique_values<TypeParam>(wt);
}



TYPED_TEST(wt_int_test, delete_)
{
    sdsl::remove(temp_file);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 3) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " test_file temp_file tmp_dir [in-memory]" << endl;
        cout << " (1) Generates a WT out of test_file; stores it in temp_file." << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    temp_file = argv[2];
    temp_dir  = argv[3];
    in_memory = argc > 4;
    if (in_memory) {
        int_vector<> data;
        load_from_file(data, test_file);
        test_file = ram_file_name(test_file);
        store_to_file(data, test_file);
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}
