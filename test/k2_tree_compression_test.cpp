        #include <sdsl/k2_tree.hpp>
#include <sdsl/k2_tree_hybrid.hpp>
#include <sdsl/k2_tree_partitioned.hpp>
#include "gtest/gtest.h"
#include "k2_tree_test_helper.hpp"

namespace {

    using namespace sdsl;
    using namespace std;

    typedef int_vector<>::size_type size_type;

    string test_file;
    string temp_file;
    bool in_memory;

    template<class T>
    class k2_tree_compression_test : public ::testing::Test {
    };

    using testing::Types;

    typedef k2_tree_hybrid<2, 2, 2, 2, bit_vector, bit_vector> hybrid_k2_2222_b_comp;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, bit_vector> hybrid_k2_4524_b_comp;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, bit_vector> hybrid_k2_4524_b;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, rrr_vector<63>> hybrid_k2_8528_b_rrr_comp;
    typedef k2_tree<2, bit_vector, bit_vector> k2;
    typedef k2_tree<4, bit_vector, bit_vector> k4;
    typedef k2_tree<4, bit_vector, bit_vector> k8;
    typedef k2_tree<4, bit_vector, bit_vector> k16;


    typedef Types<
            k8,
            k2,
            k4,
            k16,
            hybrid_k2_2222_b_comp,
            hybrid_k2_4524_b_comp,
            hybrid_k2_8528_b_rrr_comp,
            k2_tree_partitioned<k2>,
            k2_tree_partitioned<k8>,
            k2_tree_partitioned<hybrid_k2_4524_b>
    > Implementations;

    TYPED_TEST_CASE(k2_tree_compression_test, Implementations);

    TYPED_TEST(k2_tree_compression_test, CreateAndStoreTest) {
        int_vector_buffer<> buf_x(test_file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(test_file + ".y", std::ios::in);
        TypeParam k2treap(buf_x, buf_y, COUNTING_SORT);
        //construct(k2treap, test_file);
        ASSERT_TRUE(store_to_file(k2treap, temp_file));
        TypeParam k2treap2;
        ASSERT_TRUE(load_from_file(k2treap2, temp_file));
        ASSERT_EQ(k2treap, k2treap2);
    }

    /*
    TYPED_TEST(k2_tree_compression_test, direct_links) {
        TypeParam k2treap;
        perform_direct_links_test(k2treap, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.direct_links(source_id, result);
        });
    }*/



    /*
    TYPED_TEST(k2_tree_compression_test, inverse_links) {
        TypeParam k2treap;
        perform_inverse_links_test(k2treap, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.inverse_links(source_id, result);
        });
    }*/
    TYPED_TEST(k2_tree_compression_test, direct_links_2_legacy_dac) {
        TypeParam k2treap;
        perform_direct_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.direct_links2(source_id, result);
        }, LEGACY_DAC);
    }

    TYPED_TEST(k2_tree_compression_test, inverse_links2_legacy_dac) {
        TypeParam k2treap;
        perform_inverse_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.inverse_links2(source_id, result);
        }, LEGACY_DAC);
    }

    TYPED_TEST(k2_tree_compression_test, check_link_legacy_dac) {
        TypeParam k2treap;
        perform_check_link_test(k2treap, temp_file, test_file, [&k2treap](std::pair<uint64_t, uint64_t> asd) {
            return k2treap.check_link(asd);
        }, DAC);
    }

    TYPED_TEST(k2_tree_compression_test, direct_links_2_dac) {
        TypeParam k2treap;
        perform_direct_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.direct_links2(source_id, result);
        }, DAC);
    }

    TYPED_TEST(k2_tree_compression_test, inverse_links2_dac) {
        TypeParam k2treap;
        perform_inverse_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.inverse_links2(source_id, result);
        }, LEGACY_DAC);
    }

    TYPED_TEST(k2_tree_compression_test, check_link_dac) {
        TypeParam k2treap;
        perform_check_link_test(k2treap, temp_file, test_file, [&k2treap](std::pair<uint64_t, uint64_t> asd) {
            return k2treap.check_link(asd);
        }, DAC);
    }

        TYPED_TEST(k2_tree_compression_test, direct_links_2_wt_int) {
        TypeParam k2treap;
        perform_direct_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.direct_links2(source_id, result);
        }, WT_INT);
    }

    TYPED_TEST(k2_tree_compression_test, inverse_links2_wt_int) {
        TypeParam k2treap;
        perform_inverse_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.inverse_links2(source_id, result);
        }, WT_INT);
    }

    TYPED_TEST(k2_tree_compression_test, check_link_wt_int) {
        TypeParam k2treap;
        perform_check_link_test(k2treap, temp_file, test_file, [&k2treap](std::pair<uint64_t, uint64_t> asd) {
            return k2treap.check_link(asd);
        }, WT_INT);
    }
/*
    TYPED_TEST(k2_tree_compression_test, direct_links_2_wt_int_dict) {
        TypeParam k2treap;
        perform_direct_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.direct_links2(source_id, result);
        }, WT_INT_DICT);
    }

    TYPED_TEST(k2_tree_compression_test, inverse_links2_wt_int_dict) {
        TypeParam k2treap;
        perform_inverse_links_test(k2treap, temp_file, test_file, [&k2treap](uint64_t source_id, std::vector<uint64_t> &result) {
            k2treap.inverse_links2(source_id, result);
        }, WT_INT_DICT);
    }

    TYPED_TEST(k2_tree_compression_test, check_link_wt_int_dict) {
        TypeParam k2treap;
        perform_check_link_test(k2treap, temp_file, test_file, [&k2treap](std::pair<uint64_t, uint64_t> asd) {
            return k2treap.check_link(asd);
        }, WT_INT_DICT);
    }
*/
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
