#include "sdsl/k2_tree_algorithm.hpp"
#include "gtest/gtest.h"

namespace {

    using namespace sdsl;
    using namespace std;

    typedef int_vector<>::size_type size_type;

    string test_file;
    string temp_file;
    bool in_memory;

    template<class T>
    class k2_tree_access_shortcut_test : public ::testing::Test {
    };

    using testing::Types;

    typedef k2_tree_hybrid<2, 2, 2, 2, bit_vector, rrr_vector<63>, false, 2> hybrid_k2_2222_b_rrr;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, rrr_vector<63>, true, 4> hybrid_k2_4524_b_rrr;
    typedef k2_tree_hybrid<2, 5, 2, 8, bit_vector, rrr_vector<63>, false, 2> hybrid_k2_2528_b_rrr;
    typedef k2_tree_hybrid<16, 5, 2, 16, bit_vector, rrr_vector<63>, true, 2> hybrid_k2_165216_b_rrr;
    typedef k2_tree<2, bit_vector, bit_vector, true, 2> k2;
    typedef k2_tree<2, bit_vector, rrr_vector<63>, false, 4> k2rrr;

    typedef Types<
            k2,
            k2rrr,
            k2_tree_partitioned<3, k2>,
            k2_tree_partitioned<8, k2rrr>,
            hybrid_k2_2222_b_rrr,
            hybrid_k2_4524_b_rrr,
            hybrid_k2_2528_b_rrr,
            hybrid_k2_165216_b_rrr,
            k2_tree_partitioned<2,hybrid_k2_4524_b_rrr>,
            k2_tree_partitioned<3,hybrid_k2_2528_b_rrr>,
            k2_tree_partitioned<4,hybrid_k2_165216_b_rrr>
    > Implementations;

    TYPED_TEST_CASE(k2_tree_access_shortcut_test, Implementations);

    TYPED_TEST(k2_tree_access_shortcut_test, CreateAndStoreTest) {
        int_vector_buffer<> buf_x(test_file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(test_file + ".y", std::ios::in);
        try {
            TypeParam k2treap(buf_x, buf_y, false);
            //construct(k2treap, test_file);
            ASSERT_TRUE(store_to_file(k2treap, temp_file));
            TypeParam k2treap2;
            ASSERT_TRUE(load_from_file(k2treap2, temp_file));
            ASSERT_EQ(k2treap, k2treap2);

            perform_check_link_test(k2treap, [&k2treap](std::pair<uint64_t, uint64_t> asd) {
                return k2treap.check_link(asd);
            });

            perform_check_link_test(k2treap, [&k2treap](std::pair<uint64_t, uint64_t> asd) {
                return k2treap.check_link_shortcut(asd);
            });
        } catch (std::runtime_error const& e){
            std::cerr << "Exception occured " << e.what() << std::endl;
            //quite hacky comparing strings
            if (strcmp(e.what(), "shortcut size must be smaller than tree height -2") != 0){
                FAIL();
            }
        }

    }

    template<typename Function>
    void check_link_test(
            uint64_t p,
            uint64_t q,
            const int_vector<> &x,
            const int_vector<> &y,
            Function check_link
    ) {
        std::pair<uint64_t, uint64_t> asd = std::make_pair(p, q);
        bool result = check_link(asd);

        bool actual = false;
        for (uint64_t i = 0; i < x.size(); ++i) {
            if (x[i] == p) {
                if (y[i] == q) {
                    actual = true;
                }
            }
        }
        ASSERT_EQ(actual, result);
    }

    template<class t_k2treap, typename Function>
    void perform_check_link_test(t_k2treap &k2treap, Function check_link) {
        int_vector<> x, y;
        ASSERT_TRUE(load_from_file(x, test_file + ".x"));
        ASSERT_TRUE(load_from_file(y, test_file + ".y"));
        ASSERT_EQ(x.size(), y.size());
        ASSERT_EQ(x.size(), k2treap.size());
        if (x.size() > 0) {
            std::mt19937_64 rng;
            std::uniform_int_distribution<uint64_t> distribution(0, x.size() - 1);
            auto dice = bind(distribution, rng);
            //mostly negative
            for (size_t i = 0; i < 100; ++i) {
                auto idx = dice();
                uint64_t xx = x[idx];
                auto idy = dice();
                uint64_t yy = y[idy];
                check_link_test(xx, yy, x, y, check_link);
            }

            //positive
            for (size_t i = 0; i < 100; ++i) {
                auto idx = dice();
                uint64_t xx = x[idx];
                uint64_t yy = y[idx];
                check_link_test(xx, yy, x, y, check_link);
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

TYPED_TEST(k2_tree_test, count)
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

int main(int argc, char **argv) {
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
    test_file = argv[1];
    temp_file = argv[2];
    in_memory = argc > 3;
    if (in_memory) {
        auto load_and_store_in_mem = [&](string suf) {
            int_vector<> data;
            string file = temp_file + suf;
            load_vector_from_file(data, file);
            string ram_file = ram_file_name(file);
            store_to_file(data, ram_file);
        };
        load_and_store_in_mem("x");
        load_and_store_in_mem("y");
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}
