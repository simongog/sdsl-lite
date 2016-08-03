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
    class k2_tree_test : public ::testing::Test {
    };

    using testing::Types;

    typedef k2_tree_hybrid<2, 2, 2, 2, bit_vector, rrr_vector<63>> hybrid_k2_2222_b_rrr;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, rrr_vector<63>> hybrid_k2_4524_b_rrr;
    typedef k2_tree_hybrid<2, 5, 2, 8, bit_vector, rrr_vector<63>> hybrid_k2_2528_b_rrr;
    typedef k2_tree_hybrid<16, 5, 2, 16, bit_vector, rrr_vector<63>> hybrid_k2_165216_b_rrr;
    typedef k2_tree_hybrid<8, 5, 2, 3, bit_vector, rrr_vector<63>> hybrid_k2_8523_b_rrr;
    typedef k2_tree_hybrid<3, 5, 2, 4, bit_vector, rrr_vector<63>> hybrid_k2_3524_b_rrr;
    typedef k2_tree_hybrid<2, 2, 2, 2, bit_vector, bit_vector, true> hybrid_k2_2222_b_comp;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, bit_vector, true> hybrid_k2_4524_b_comp;
    typedef k2_tree_hybrid<16, 5, 2, 16, bit_vector, rrr_vector<63>, true> hybrid_k2_165216_b_rrr_comp;
    typedef k2_tree<2, bit_vector> k2;
    typedef k2_tree<2, rrr_vector<63>> k2rrr;
    typedef k2_tree<3, bit_vector> k3;
    typedef k2_tree<4, bit_vector> k4;
    typedef k2_tree<6, bit_vector> k6;
    typedef k2_tree<8, bit_vector> k8;
    typedef k2_tree<16, bit_vector> k16;
    typedef k2_tree<2, bit_vector, bit_vector, true> k2comp;
    typedef k2_tree<4, bit_vector, bit_vector, true> k4comp;
    typedef k2_tree<8, bit_vector, bit_vector, true> k8comp;
    typedef k2_tree<16, bit_vector, bit_vector, true> k16comp;


    typedef Types<
            k8comp/*,
            k2comp,
            k4comp,
            k16comp,
            k2,
            k2rrr,
            k3,
            k4,
            k6,
            k8,
            k16,
            hybrid_k2_2222_b_comp,
            hybrid_k2_4524_b_comp,
            hybrid_k2_165216_b_rrr_comp,
            k2_tree_partitioned<2, k2>,
            k2_tree_partitioned<4, k2rrr>,
            k2_tree_partitioned<3, k3>,
            k2_tree_partitioned<8, k4>,
            k2_tree_partitioned<16,k6>,
            k2_tree_partitioned<16,k8>,
            k2_tree_partitioned<16,k16>,
            k2_tree_partitioned<2, k2comp>,
            k2_tree_partitioned<4, k8comp>,
            hybrid_k2_2222_b_rrr,
            hybrid_k2_4524_b_rrr,
            hybrid_k2_2528_b_rrr,
            hybrid_k2_165216_b_rrr,
            hybrid_k2_8523_b_rrr,
            hybrid_k2_3524_b_rrr,
            k2_tree_partitioned<2,hybrid_k2_4524_b_rrr>,
            k2_tree_partitioned<3,hybrid_k2_2528_b_rrr>,
            k2_tree_partitioned<4,hybrid_k2_165216_b_rrr>,
            k2_tree_partitioned<8,hybrid_k2_8523_b_rrr>,
            k2_tree_partitioned<16,hybrid_k2_3524_b_rrr>*/
    > Implementations;

    TYPED_TEST_CASE(k2_tree_test, Implementations);

    TYPED_TEST(k2_tree_test, CreateAndStoreTest) {
        TypeParam k2treap;
        construct(k2treap, test_file);
        ASSERT_TRUE(store_to_file(k2treap, temp_file));
    }

    TYPED_TEST(k2_tree_test, ConstructCompareTest) {
        TypeParam k2treap;
        TypeParam k2treap2;
        construct(k2treap, test_file);
        construct_bottom_up(k2treap2, test_file);
        std::cout << "Comparing Results" << std::endl;

        if (!( k2treap == k2treap2)){
            std::cout << "Results differ" << std::endl;
        }

        ASSERT_EQ(k2treap, k2treap2);
    }

    TYPED_TEST(k2_tree_test, size) {
        TypeParam k2treap;
        ASSERT_TRUE(load_from_file(k2treap, temp_file));
        int_vector<> x, y;
        ASSERT_TRUE(load_from_file(x, test_file + ".x"));
        ASSERT_TRUE(load_from_file(y, test_file + ".y"));
        ASSERT_EQ(x.size(), y.size());
        ASSERT_EQ(x.size(), k2treap.size());
    }

    template<typename t_k2treap, typename Function>
    void direct_links_test(
            uint64_t source_id,
            const int_vector<> &x,
            const int_vector<> &y,
            Function direct_links,
            t_k2treap &k2treap
    ) {
        std::vector<uint64_t> result;
        direct_links(source_id, result);

        typedef tuple<uint64_t, uint64_t> t_xy;
        vector<t_xy> vec;
        for (uint64_t i = 0; i < x.size(); ++i) {
            if (x[i] == source_id) {
                vec.emplace_back(x[i], y[i]);
            }
        }
        sort(vec.begin(), vec.end(), [](const t_xy &a, const t_xy &b) {
            if (get<0>(a) != get<0>(b))
                return get<0>(a) < get<0>(b);
            return get<1>(a) < get<1>(b);
        });


        if (result.size() != vec.size()){
            std::cout << "Source_id: " << source_id << std::endl;
            std::cout << "Result:" << std::endl;
            for (auto asd : result){
                std::cout << asd << std::endl;
            }

            std::cout << "Vec:" << std::endl;
            for (auto asd : vec){
                std::cout << get<1>(asd) << std::endl;
            }
/*
            std::cout << "X, Y vectors" << std::endl;
            for (uint64_t i = 0; i < x.size(); ++i) {
                std::cout << x[i] << "," << y[i] << std::endl;
            }
*/
            std::vector<uint64_t> result;
            k2treap.direct_links2(source_id, result);
        }

        ASSERT_EQ(result.size(), vec.size());

        uint64_t cnt = 0;
        std::vector<uint64_t>::iterator res_it = result.begin();
        for (auto direct_link: result) {
            ASSERT_TRUE(cnt < vec.size());
            ASSERT_EQ(get<1>(vec[cnt]), direct_link);
            ++res_it;
            ++cnt;
        }
    }

    template<class t_k2treap, typename Function>
    void perform_direct_links_test(t_k2treap &k2treap, Function direct_links) {

        ASSERT_TRUE(load_from_file(k2treap, temp_file));
        int_vector<> x, y;
        ASSERT_TRUE(load_from_file(x, test_file + ".x"));
        ASSERT_TRUE(load_from_file(y, test_file + ".y"));
        ASSERT_EQ(x.size(), y.size());
        ASSERT_EQ(x.size(), k2treap.size());
        if (x.size() > 0) {
            std::mt19937_64 rng;
            std::uniform_int_distribution<uint64_t> distribution(0, x.size() - 1);
            auto dice = bind(distribution, rng);
            for (size_t i = 0; i < 100; ++i) {
                auto idx = dice();
                uint64_t xx = x[idx];
                direct_links_test(xx, x, y, direct_links, k2treap);
            }
        }
    }

    /*
    TYPED_TEST(k2_tree_test, direct_links) {
        TypeParam k2treap;
        perform_direct_links_test(k2treap, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.direct_links(source_id, result);
        });
    }*/

    TYPED_TEST(k2_tree_test, direct_links_2) {
        TypeParam k2treap;
        perform_direct_links_test(k2treap, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.direct_links2(source_id, result);
        });
    }

    template<typename Function>
    void inverse_links_test(
            uint64_t source_id,
            const int_vector<> &x,
            const int_vector<> &y,
            Function inverse_links
    ) {
        std::vector<uint64_t> result;
        inverse_links(source_id, result);

        typedef tuple<uint64_t, uint64_t> t_xy;
        vector<t_xy> vec;
        for (uint64_t i = 0; i < x.size(); ++i) {
            if (y[i] == source_id) {
                vec.emplace_back(x[i], y[i]);
            }
        }
        sort(vec.begin(), vec.end(), [](const t_xy &a, const t_xy &b) {
            if (get<0>(a) != get<0>(b))
                return get<0>(a) < get<0>(b);
            return get<1>(a) < get<1>(b);
        });

        ASSERT_EQ(result.size(), vec.size());

        uint64_t cnt = 0;
        std::vector<uint64_t>::iterator res_it = result.begin();
        for (auto link_source: result) {
            ASSERT_TRUE(cnt < vec.size());
            ASSERT_EQ(get<0>(vec[cnt]), link_source);
            ++res_it;
            ++cnt;
        }
    }

    template<class t_k2treap, typename Function>
    void perform_inverse_links_test(t_k2treap &k2treap, Function inverse_links) {

        ASSERT_TRUE(load_from_file(k2treap, temp_file));
        int_vector<> x, y;
        ASSERT_TRUE(load_from_file(x, test_file + ".x"));
        ASSERT_TRUE(load_from_file(y, test_file + ".y"));
        ASSERT_EQ(x.size(), y.size());
        ASSERT_EQ(x.size(), k2treap.size());
        if (x.size() > 0) {
            std::mt19937_64 rng;
            std::uniform_int_distribution<uint64_t> distribution(0, x.size() - 1);
            auto dice = bind(distribution, rng);
            for (size_t i = 0; i < 100; ++i) {
                auto idx = dice();
                uint64_t xx = x[idx];
                inverse_links_test(xx, x, y, inverse_links);
            }
        }
    }

    /*
    TYPED_TEST(k2_tree_test, inverse_links) {
        TypeParam k2treap;
        perform_inverse_links_test(k2treap, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.inverse_links(source_id, result);
        });
    }*/

    TYPED_TEST(k2_tree_test, inverse_links2) {
        TypeParam k2treap;
        perform_inverse_links_test(k2treap, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.inverse_links2(source_id, result);
        });
    }

    template<class t_k2treap>
    void check_link_test(
            t_k2treap &k2treap,
            uint64_t p,
            uint64_t q,
            const int_vector<> &x,
            const int_vector<> &y
    ) {
        std::pair<uint64_t, uint64_t> asd = std::make_pair(p, q);
        bool result = k2treap.check_link(asd);

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

    /*
    TYPED_TEST(k2_tree_test, check_link) {
        TypeParam k2treap;
        ASSERT_TRUE(load_from_file(k2treap, temp_file));
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
            for (size_t i = 0; i < 50; ++i) {
                auto idx = dice();
                uint64_t xx = x[idx];
                auto idy = dice();
                uint64_t yy = y[idy];
                check_link_test(k2treap, xx, yy, x, y);
            }

            //postive
            for (size_t i = 0; i < 50; ++i) {
                auto idx = dice();
                uint64_t xx = x[idx];
                uint64_t yy = y[idx];
                check_link_test(k2treap, xx, yy, x, y);
            }
        }
    }*/

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
        load_and_store_in_mem("w");
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}
