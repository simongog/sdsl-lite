#include "sdsl/wt_topk.hpp"
#include "sdsl/wt_topk_algorithm.hpp"
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
typedef tuple<uint64_t, uint64_t, uint64_t> t_xyw;

string test_file;
string temp_file;
bool in_memory;

template<class T>
class WtTopkTest : public ::testing::Test { };

using testing::Types;

typedef Types<
<<<<<<< HEAD
wt_topk<wt_int<>, rmq_succinct_sct<false>, dac_vector<> >
=======
wt_topk<>
>>>>>>> wt_topk
> Implementations;

TYPED_TEST_CASE(WtTopkTest, Implementations);

TYPED_TEST(WtTopkTest, CreateAndStoreTest)
{
<<<<<<< HEAD
    TypeParam topk_wt;
    construct(topk_wt, test_file);
    ASSERT_TRUE(store_to_file(topk_wt, temp_file));
}

struct my_xyw_comp {
    bool operator()(const t_xyw& a, const t_xyw& b) const
    {
        if (get<2>(a) != get<2>(b))
            return get<2>(a) > get<2>(b);
        else if (get<0>(a) != get<0>(b))
            return get<0>(a) < get<0>(b);
        return get<1>(a) < get<1>(b);
    }
};

template<class t_topk_wt>
void topk_test(
    const t_topk_wt& topk_wt,
=======
    TypeParam wttopk;
    construct(wttopk, test_file);
    ASSERT_TRUE(store_to_file(wttopk, temp_file));
}


template<class t_wttopk>
void topk_test(
    const t_wttopk& wttopk,
>>>>>>> wt_topk
    complex<uint64_t> min_xy,
    complex<uint64_t> max_xy,
    const int_vector<>& x,
    const int_vector<>& y,
    const int_vector<>& w)
{
<<<<<<< HEAD
    auto res_it = top_k(topk_wt, {real(min_xy),imag(min_xy)}, {real(max_xy),imag(max_xy)});
=======
    auto res_it = top_k(wttopk, {real(min_xy),imag(min_xy)}, {real(max_xy),imag(max_xy)});
    typedef tuple<uint64_t, uint64_t, uint64_t> t_xyw;
>>>>>>> wt_topk
    vector<t_xyw> vec;
    for (uint64_t i = 0; i < x.size(); ++i) {
        if (x[i] >= real(min_xy) and x[i] <= real(max_xy)
            and y[i] >= imag(min_xy) and y[i] <= imag(max_xy)) {
            vec.emplace_back(x[i], y[i], w[i]);
        }
    }
<<<<<<< HEAD
    sort(vec.begin(), vec.end(), my_xyw_comp());
    uint64_t cnt = 0;
    vector<t_xyw> vec2;
    while (res_it) {
        ASSERT_TRUE(cnt < vec.size());
        auto p = *res_it;
        vec2.emplace_back(real(p.first), imag(p.first), p.second);
        if (vec2.size() > 1) {
            EXPECT_TRUE(get<2>(vec2[vec2.size()-2]) >= get<2>(vec2[vec2.size()-1]))
                    << get<2>(vec2[vec2.size()-2]) <<" < " << get<2>(vec2[vec2.size()-1]);
        }
        ++res_it;
        ++cnt;
=======
    sort(vec.begin(), vec.end(), [](const t_xyw& a, const t_xyw& b) {
            return get<2>(a) > get<2>(b);
    });

    uint64_t cnt = 0;
    if (vec.size() != 0) {
        while (res_it) {
            ASSERT_TRUE(cnt < vec.size());
            auto p = *res_it;
            ASSERT_EQ(get<2>(vec[cnt]), p.second);
            ++res_it;
            ++cnt;
        }
>>>>>>> wt_topk
    }
    ASSERT_FALSE(res_it);
    sort(vec2.begin(), vec2.end(), my_xyw_comp());
    ASSERT_EQ(vec.size(), vec2.size());
    for (size_t i=0; i<vec.size(); ++i) {
        EXPECT_EQ(vec[i], vec2[i]) << "i="<<i<<" min_xy="<<min_xy<<" max_xy="<<max_xy;
    }
}

TYPED_TEST(WtTopkTest, SizeAndTopk)
{
<<<<<<< HEAD
    TypeParam topk_wt;
    ASSERT_TRUE(load_from_file(topk_wt, temp_file));
=======
    TypeParam wttopk;
    ASSERT_TRUE(load_from_file(wttopk, temp_file));
>>>>>>> wt_topk
    int_vector<> x,y,w;
    ASSERT_TRUE(load_from_file(x, test_file+".x"));
    ASSERT_TRUE(load_from_file(y, test_file+".y"));
    ASSERT_EQ(x.size(), y.size());
    ASSERT_TRUE(load_from_file(w, test_file+".w"));
    ASSERT_EQ(x.size(), w.size());
<<<<<<< HEAD
    ASSERT_EQ(x.size(), topk_wt.size());
=======
    ASSERT_EQ(x.size(), wttopk.size());
>>>>>>> wt_topk
    uint64_t maxx=0, maxy=0;
    if (x.size() > 0) {
        maxx =  *max_element(x.begin(), x.end());
        maxy =  *max_element(y.begin(), y.end());
    }
    uint64_t minx=0, miny=0;
<<<<<<< HEAD
    topk_test(topk_wt, {minx,miny}, {maxx,maxy}, x, y, w);
=======
    topk_test(wttopk, {minx,maxx}, {miny,maxy}, x, y, w);
>>>>>>> wt_topk

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
            if (y >= dd)
                miny = yy - dd;
<<<<<<< HEAD
            topk_test(topk_wt, {minx, miny}, {maxx,maxy}, x, y, w);
=======
            topk_test(wttopk, {minx, miny}, {maxx,maxy}, x, y, w);
>>>>>>> wt_topk
        }
    }
}
// RK: Commenting this one, as we may implement range_3d
//
//template<class t_wttopk>
//void range3d_test(
//    const t_wttopk& wttopk,
//    complex<uint64_t> min_xy,
//    complex<uint64_t> max_xy,
//    complex<uint64_t> z,
//    const int_vector<>& x,
//    const int_vector<>& y,
//    const int_vector<>& w)
//{
//    auto res_it = range_3d(wttopk, {real(min_xy),imag(min_xy)},
//    {real(max_xy),imag(max_xy)},
//    {real(z), imag(z)});
//    typedef tuple<uint64_t, uint64_t, uint64_t> t_xyw;
//    vector<t_xyw> vec;
//    for (uint64_t i = 0; i < x.size(); ++i) {
//        if (x[i] >= real(min_xy) and x[i] <= real(max_xy)
//            and y[i] >= imag(min_xy) and y[i] <= imag(max_xy)) {
//            vec.emplace_back(x[i], y[i], w[i]);
//        }
//    }
//    sort(vec.begin(), vec.end(), [](const t_xyw& a, const t_xyw& b) {
//        if (get<2>(a) != get<2>(b))
//            return get<2>(a) > get<2>(b);
//        else if (get<0>(a) != get<0>(b))
//            return get<0>(a) < get<0>(b);
//        return get<1>(a) < get<1>(b);
//    });
//    uint64_t cnt = 0;
//    while (res_it) {
//        ASSERT_TRUE(cnt < vec.size());
//        auto p = *res_it;
//        ASSERT_EQ(get<2>(vec[cnt]), p.second);
//        ASSERT_EQ(get<0>(vec[cnt]), real(p.first));
//        ASSERT_EQ(get<1>(vec[cnt]), imag(p.first));
//        ++res_it;
//        ++cnt;
//    }
//    ASSERT_FALSE(res_it);
//}
//
//TYPED_TEST(WtTopkTest, Range3d)
//{
//    TypeParam wttopk;
//    ASSERT_TRUE(load_from_file(wttopk, temp_file));
//    int_vector<> x,y,w;
//    ASSERT_TRUE(load_from_file(x, test_file+".x"));
//    ASSERT_TRUE(load_from_file(y, test_file+".y"));
//    ASSERT_EQ(x.size(), y.size());
//    ASSERT_TRUE(load_from_file(w, test_file+".w"));
//    ASSERT_EQ(x.size(), w.size());
//    ASSERT_EQ(x.size(), k2treap.size());
//    if (x.size() > 0) {
//        std::mt19937_64 rng;
//        std::uniform_int_distribution<uint64_t> distribution(0, x.size()-1);
//        auto dice = bind(distribution, rng);
//        for (size_t i=0; i<20; ++i) {
//            auto idx = dice();
//            uint64_t xx = x[idx];
//            uint64_t yy = y[idx];
//            uint64_t ww = w[idx];
//            uint64_t dd = 20;
//            uint64_t dw = 100;
//            uint64_t minx=0, miny=0, maxx=xx+dd, maxy=yy+dd, minw=0, maxw=ww+dw;
//            if (xx >= dd)
//                minx = xx - dd;
//            if (yy >= dd)
//                miny = yy - dd;
//            if (ww >= dw)
//                minw = ww - dw;
//            range3d_test(wttopk, {minx, miny}, {maxx,maxy}, {minw,maxw}, x, y, w);
//        }
//    }
//}

<<<<<<< HEAD
/*

template<class t_topk_wt>
void range3d_test(
    const t_topk_wt& topk_wt,
    complex<uint64_t> min_xy,
    complex<uint64_t> max_xy,
    complex<uint64_t> z,
    const int_vector<>& x,
    const int_vector<>& y,
    const int_vector<>& w)
{
    auto res_it = range_3d(topk_wt, {real(min_xy),imag(min_xy)},
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
*/

/*
TYPED_TEST(WtTopkTest, Range3d)
{
    TypeParam topk_wt;
    ASSERT_TRUE(load_from_file(topk_wt, temp_file));
    int_vector<> x,y,w;
    ASSERT_TRUE(load_from_file(x, test_file+".x"));
    ASSERT_TRUE(load_from_file(y, test_file+".y"));
    ASSERT_EQ(x.size(), y.size());
    ASSERT_TRUE(load_from_file(w, test_file+".w"));
    ASSERT_EQ(x.size(), w.size());
    ASSERT_EQ(x.size(), topk_wt.size());
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
            range3d_test(topk_wt, {minx, miny}, {maxx,maxy}, {minw,maxw}, x, y, w);
        }
    }
}

*/

template<class t_topk_wt>
void count_test(
    const t_topk_wt& topk_wt,
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
    ASSERT_EQ(cnt, count(topk_wt, {real(min_xy),imag(min_xy)}, {real(max_xy),imag(max_xy)}));
}

TYPED_TEST(WtTopkTest, Count)
{
    TypeParam topk_wt;
    ASSERT_TRUE(load_from_file(topk_wt, temp_file));
    int_vector<> x,y;
    ASSERT_TRUE(load_from_file(x, test_file+".x"));
    ASSERT_TRUE(load_from_file(y, test_file+".y"));
    ASSERT_EQ(x.size(), y.size());
    ASSERT_EQ(x.size(), topk_wt.size());
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
            count_test(topk_wt, {std::min(x1,x2), std::min(y1,y2)}, {std::max(x1,x2),std::max(y1,y2)}, x, y);
        }
    }
}
=======
>>>>>>> wt_topk

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
