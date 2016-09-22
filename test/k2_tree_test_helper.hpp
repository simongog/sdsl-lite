//
// Created by d056848 on 8/11/16.
//

#ifndef INCLUDED_SDSL_K2_TREE_TEST_HELPER
#define INCLUDED_SDSL_K2_TREE_TEST_HELPER

#include "gtest/gtest.h"
#include "sdsl/k2_tree_algorithm.hpp"
#include "sdsl/vectors.hpp"
#include "sdsl/io.hpp"

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;


template<typename Function>
void direct_links_test(
        uint64_t source_id,
        const int_vector<> &x,
        const int_vector<> &y,
        Function direct_links
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


    if (result.size() != vec.size()) {
        std::cout << "Source_id: " << source_id << std::endl;
        std::cout << "Result:" << std::endl;
        for (auto asd : result) {
            std::cout << asd << std::endl;
        }

        std::cout << "Vec:" << std::endl;
        for (auto asd : vec) {
            std::cout << get<1>(asd) << std::endl;
        }
/*
            std::cout << "X, Y vectors" << std::endl;
            for (uint64_t i = 0; i < x.size(); ++i) {
                std::cout << x[i] << "," << y[i] << std::endl;
            }
*/
        std::vector<uint64_t> result;
        direct_links(source_id, result);
        std::cout << "Result after second try:" << std::endl;
        for (auto asd : result) {
            std::cout << asd << std::endl;
        }
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
void perform_direct_links_test(t_k2treap &k2treap, std::string temp_file, std::string test_file, Function direct_links, leaf_compression_type compression=UNCOMPRESSED) {

    ASSERT_TRUE(load_from_file(k2treap, temp_file));
    int_vector<> x, y;

    if (compression != UNCOMPRESSED){
        k2treap.compress_leaves(compression);
    }

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
            direct_links_test(xx, x, y, direct_links);
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

    if (actual != result){
        std::cout << "Result differs for: " << p <<"," << q << std::endl;
        std::cout << "Actual: " << actual << std::endl;
        std::cout << "Result: " << result << std::endl;
        check_link(asd);
    }

    ASSERT_EQ(actual, result);
}

template<class t_k2treap, typename Function>
void perform_check_link_test(t_k2treap &k2treap, std::string temp_file, std::string test_file, Function check_link, leaf_compression_type compression=UNCOMPRESSED) {
    ASSERT_TRUE(load_from_file(k2treap, temp_file));
    if (compression != UNCOMPRESSED){
        k2treap.compress_leaves(compression);
    }
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
void perform_inverse_links_test(t_k2treap &k2treap, std::string temp_file, std::string test_file, Function inverse_links, leaf_compression_type compression=UNCOMPRESSED) {
    ASSERT_TRUE(load_from_file(k2treap, temp_file));
    if (compression != UNCOMPRESSED){
        k2treap.compress_leaves(compression);
    }
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


#endif

