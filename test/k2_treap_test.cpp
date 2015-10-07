#include "sdsl/k2_treap.hpp"
#include "sdsl/bit_vectors.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <tuple>
#include <string>
#include <algorithm> // for std::min. std::sort
#include <random>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;

string test_file;
string temp_file;
bool in_memory;

template<class T>
class k2_treap_test : public ::testing::Test { };

using testing::Types;

typedef Types<
k2_treap<2, bit_vector>,
         k2_treap<2, rrr_vector<63>>,
         k2_treap<3, bit_vector>,
         k2_treap<4, rrr_vector<63>>,
         k2_treap<5, rrr_vector<63>>,
         k2_treap<6, rrr_vector<63>>,
         k2_treap<16, rrr_vector<63>>
         > Implementations;

TYPED_TEST_CASE(k2_treap_test, Implementations);

TYPED_TEST(k2_treap_test, CreateAndStoreTest)
{
    TypeParam k2treap;
    construct(k2treap, test_file);
    ASSERT_TRUE(store_to_file(k2treap, temp_file));
}

template<class t_k2treap>
void topk_test(
    const t_k2treap& k2treap,
    complex<uint64_t> min_xy,
    complex<uint64_t> max_xy,
    const int_vector<>& x,
    const int_vector<>& y,
    const int_vector<>& w)
{
    auto res_it = top_k(k2treap, {real(min_xy),imag(min_xy)}, {real(max_xy),imag(max_xy)});
    typedef tuple<uint64_t, uint64_t, uint64_t> t_xyw;
    vector<t_xyw> vec;
    for (uint64_t i = 0; i < x.size(); ++i) {
        if (x[i] >= real(min_xy) and x[i] <= real(max_xy)
            and y[i] >= imag(min_xy) and y[i] <= imag(max_xy)) {
            vec.emplace_back(x[i], y[i], w[i]);
        }
    }
    sort(vec.begin(), vec.end(), [](const t_xyw& a, const t_xyw& b) {
        if (get<2>(a) != get<2>(b))
            return get<2>(a) > get<2>(b);
        else if (get<0>(a) != get<0>(b))
            return get<0>(a) < get<0>(b);
        return get<1>(a) < get<1>(b);
    });
    uint64_t cnt = 0;
    while (res_it) {
        ASSERT_TRUE(cnt < vec.size());
        auto p = *res_it;
        ASSERT_EQ(get<2>(vec[cnt]), p.second);
        ASSERT_EQ(get<0>(vec[cnt]), real(p.first));
        ASSERT_EQ(get<1>(vec[cnt]), imag(p.first));
        ++res_it;
        ++cnt;
    }
    ASSERT_FALSE(res_it);
}

TYPED_TEST(k2_treap_test, size_and_top_k)
{
    TypeParam k2treap;
    ASSERT_TRUE(load_from_file(k2treap, temp_file));
    int_vector<> x,y,w;
    ASSERT_TRUE(load_from_file(x, test_file+".x"));
    ASSERT_TRUE(load_from_file(y, test_file+".y"));
    ASSERT_EQ(x.size(), y.size());
    ASSERT_TRUE(load_from_file(w, test_file+".w"));
    ASSERT_EQ(x.size(), w.size());
    ASSERT_EQ(x.size(), k2treap.size());
    uint64_t maxx=0, maxy=0;
    if (x.size() > 0) {
        maxx =  *max_element(x.begin(), x.end());
        maxy =  *max_element(y.begin(), y.end());
    }
    uint64_t minx=0, miny=0;
    topk_test(k2treap, {minx,maxx}, {miny,maxy}, x, y, w);

    if (x.size() > 0) {
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, x.size()-1);
        auto dice = bind(distribution, rng);
        for (size_t i=0; i<20; ++i) {
            auto idx = dice();
            uint64_t xx = x[idx];
            uint64_t yy = y[idx];
            uint64_t dd = 20;
            uint64_t minx=0, miny=0, maxx=xx+dd, maxy=yy+dd;
            if (xx >= dd)
                minx = xx - dd;
            if (yy >= dd)
                miny = yy - dd;
            topk_test(k2treap, {minx, miny}, {maxx,maxy}, x, y, w);
        }
    }
}

template<class t_k2treap>
void range3d_test(
    const t_k2treap& k2treap,
    complex<uint64_t> min_xy,
    complex<uint64_t> max_xy,
    complex<uint64_t> z,
    const int_vector<>& x,
    const int_vector<>& y,
    const int_vector<>& w)
{
    auto res_it = range_3d(k2treap, {real(min_xy),imag(min_xy)},
    {real(max_xy),imag(max_xy)},
    {real(z), imag(z)});
    typedef tuple<uint64_t, uint64_t, uint64_t> t_xyw;
    vector<t_xyw> vec;
    for (uint64_t i = 0; i < x.size(); ++i) {
        if (x[i] >= real(min_xy) and x[i] <= real(max_xy)
            and y[i] >= imag(min_xy) and y[i] <= imag(max_xy)) {
            vec.emplace_back(x[i], y[i], w[i]);
        }
    }
    sort(vec.begin(), vec.end(), [](const t_xyw& a, const t_xyw& b) {
        if (get<2>(a) != get<2>(b))
            return get<2>(a) > get<2>(b);
        else if (get<0>(a) != get<0>(b))
            return get<0>(a) < get<0>(b);
        return get<1>(a) < get<1>(b);
    });
    uint64_t cnt = 0;
    while (res_it) {
        ASSERT_TRUE(cnt < vec.size());
        auto p = *res_it;
        ASSERT_EQ(get<2>(vec[cnt]), p.second);
        ASSERT_EQ(get<0>(vec[cnt]), real(p.first));
        ASSERT_EQ(get<1>(vec[cnt]), imag(p.first));
        ++res_it;
        ++cnt;
    }
    ASSERT_FALSE(res_it);
}

TYPED_TEST(k2_treap_test, range_3d)
{
    TypeParam k2treap;
    ASSERT_TRUE(load_from_file(k2treap, temp_file));
    int_vector<> x,y,w;
    ASSERT_TRUE(load_from_file(x, test_file+".x"));
    ASSERT_TRUE(load_from_file(y, test_file+".y"));
    ASSERT_EQ(x.size(), y.size());
    ASSERT_TRUE(load_from_file(w, test_file+".w"));
    ASSERT_EQ(x.size(), w.size());
    ASSERT_EQ(x.size(), k2treap.size());
    if (x.size() > 0) {
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, x.size()-1);
        auto dice = bind(distribution, rng);
        for (size_t i=0; i<20; ++i) {
            auto idx = dice();
            uint64_t xx = x[idx];
            uint64_t yy = y[idx];
            uint64_t ww = w[idx];
            uint64_t dd = 20;
            uint64_t dw = 100;
            uint64_t minx=0, miny=0, maxx=xx+dd, maxy=yy+dd, minw=0, maxw=ww+dw;
            if (xx >= dd)
                minx = xx - dd;
            if (yy >= dd)
                miny = yy - dd;
            if (ww >= dw)
                minw = ww - dw;
            range3d_test(k2treap, {minx, miny}, {maxx,maxy}, {minw,maxw}, x, y, w);
        }
    }
}

template<class t_k2treap>
void count_test(
    const t_k2treap& k2treap,
    complex<uint64_t> min_xy,
    complex<uint64_t> max_xy,
    const int_vector<>& x,
    const int_vector<>& y)
{
    uint64_t cnt = 0;
    for (uint64_t i = 0; i < x.size(); ++i) {
        if (x[i] >= real(min_xy) and x[i] <= real(max_xy)
            and y[i] >= imag(min_xy) and y[i] <= imag(max_xy)) {
            ++cnt;
        }
    }
    ASSERT_EQ(cnt, count(k2treap, {real(min_xy),imag(min_xy)}, {real(max_xy),imag(max_xy)}));
}

TYPED_TEST(k2_treap_test, count)
{
    TypeParam k2treap;
    ASSERT_TRUE(load_from_file(k2treap, temp_file));
    int_vector<> x,y;
    ASSERT_TRUE(load_from_file(x, test_file+".x"));
    ASSERT_TRUE(load_from_file(y, test_file+".y"));
    ASSERT_EQ(x.size(), y.size());
    ASSERT_EQ(x.size(), k2treap.size());
    if (x.size() > 0) {
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, x.size()-1);
        auto dice = bind(distribution, rng);
        for (size_t i=0; i<3; ++i) {
            auto idx1 = dice();
            auto idx2 = dice();
            uint64_t x1 = x[idx1];
            uint64_t y1 = y[idx1];
            uint64_t x2 = x[idx2];
            uint64_t y2 = y[idx2];
            count_test(k2treap, {std::min(x1,x2), std::min(y1,y2)}, {std::max(x1,x2),std::max(y1,y2)}, x, y);
        }
    }
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 3) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " file temp_file [in-memory]" << endl;
        cout << " (1) Generates a k2-treap out of file.x, file.y, and file.w." << endl;
        cout << "     Result is stored in temp_file." << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file    = argv[1];
    temp_file    = argv[2];
    in_memory    = argc > 3;
    if (in_memory) {
        auto load_and_store_in_mem = [&](string suf) {
            int_vector<> data;
            string file = temp_file + suf;
            load_vector_from_file(data,file);
            string ram_file = ram_file_name(file);
            store_to_file(data, ram_file);
        };
        load_and_store_in_mem("x");
        load_and_store_in_mem("y");
        load_and_store_in_mem("w");
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}
