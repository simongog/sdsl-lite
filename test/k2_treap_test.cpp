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

TYPED_TEST(k2_treap_test, ConstructCompareTest)
{
        TypeParam k2treap;
        TypeParam k2treap2;
        construct(k2treap, test_file);
        construct_bottom_up(k2treap2, test_file);
        ASSERT_EQ(k2treap,k2treap2);
}

TYPED_TEST(k2_treap_test, size)
{
    TypeParam k2treap;
    ASSERT_TRUE(load_from_file(k2treap, temp_file));
    int_vector<> x,y;
    ASSERT_TRUE(load_from_file(x, test_file+".x"));
    ASSERT_TRUE(load_from_file(y, test_file+".y"));
    ASSERT_EQ(x.size(), y.size());
    ASSERT_EQ(x.size(), k2treap.size());
}

template<class t_k2treap>
void direct_links_test(
    const t_k2treap& k2treap,
    uint64_t source_id,
    const int_vector<>& x,
    const int_vector<>& y
    )
{
    std::vector<uint64_t> result;
    k2treap.direct_links(source_id, result);

    typedef tuple<uint64_t, uint64_t> t_xy;
    vector<t_xy> vec;
    for (uint64_t i = 0; i < x.size(); ++i) {
        if (x[i] == source_id) {
            vec.emplace_back(x[i], y[i]);
        }
    }
    sort(vec.begin(), vec.end(), [](const t_xy& a, const t_xy& b) {
        if (get<0>(a) != get<0>(b))
            return get<0>(a) < get<0>(b);
        return get<1>(a) < get<1>(b);
    });

    uint64_t cnt = 0;
    std::vector<uint64_t>::iterator res_it = result.begin();
    for (auto direct_link: result){
        ASSERT_TRUE(cnt < vec.size());
        ASSERT_EQ(get<1>(vec[cnt]), direct_link);
        ++res_it;
        ++cnt;
    }
}

TYPED_TEST(k2_treap_test, direct_links)
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
        for (size_t i=0; i<20; ++i) {
            auto idx = dice();
            uint64_t xx = x[idx];
            direct_links_test(k2treap, xx, x, y);
        }
    }
}

    /*
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
*/

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
