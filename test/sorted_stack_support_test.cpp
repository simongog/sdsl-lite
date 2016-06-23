#include "sdsl/sorted_stack_support.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <string>
#include <random>

namespace
{

std::string temp_dir;

// The fixture for testing class sorted_stack_support.
class sorted_stack_support_test : public ::testing::Test
{
    protected:

        sorted_stack_support_test() {}

        virtual ~sorted_stack_support_test() {}

        virtual void SetUp() {}

        virtual void TearDown() {}
};

void compare_stacks(std::stack<uint64_t>& exp, sdsl::sorted_stack_support& sss)
{
    ASSERT_EQ(exp.size(), sss.size());
    std::stack<uint64_t> tmp;
    while (!exp.empty()) {
        ASSERT_EQ(exp.top(), sss.top());
        tmp.push(exp.top());
        exp.pop();
        sss.pop();
    }
    ASSERT_EQ(exp.size(), sss.size());
    // Restore stacks
    while (!tmp.empty()) {
        exp.push(tmp.top());
        sss.push(tmp.top());
        tmp.pop();
    }
}

//! Test Constructors
TEST_F(sorted_stack_support_test, constructors)
{
    static_assert(std::is_copy_constructible<sdsl::sorted_stack_support>::value, "Type is not copy constructible");
    static_assert(std::is_move_constructible<sdsl::sorted_stack_support>::value, "Type is not move constructible");
    static_assert(std::is_copy_assignable<sdsl::sorted_stack_support>::value, "Type is not copy assignable");
    static_assert(std::is_move_assignable<sdsl::sorted_stack_support>::value, "Type is not move assignable");
    std::stack<uint64_t> exp;
    sdsl::sorted_stack_support sss1(100000+10);
    {
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, 10);
        auto dice = bind(distribution, rng);
        for (uint64_t k=0; k<100000; ++k) {
            uint64_t value = k+dice();
            if (exp.empty() or exp.top() < value) {
                exp.push(value);
                sss1.push(value);
            }
        }
    }

    // Copy-constructor
    sdsl::sorted_stack_support sss2(sss1);
    compare_stacks(exp, sss2);

    // Move-constructor
    sdsl::sorted_stack_support sss3(std::move(sss2));
    compare_stacks(exp, sss3);
    // ASSERT_EQ((uint64_t)0, sss2.size());

    // Copy-Assign
    sdsl::sorted_stack_support sss4(0);
    sss4 = sss3;
    compare_stacks(exp, sss4);

    // Move-Assign
    sdsl::sorted_stack_support sss5(0);
    sss5 = std::move(sss4);
    compare_stacks(exp, sss5);
    // ASSERT_EQ((uint64_t)0, sss4.size());
}

//! Test Operations push, top and pop
TEST_F(sorted_stack_support_test, push_top_and_pop)
{
    for (uint64_t i=0; i<20; ++i) {
        std::mt19937_64 rng(i);
        std::uniform_int_distribution<uint64_t> distribution(0, i*i);
        auto dice = bind(distribution, rng);
        std::stack<uint64_t> exp;
        sdsl::sorted_stack_support sss(1000000+i*i);
        ASSERT_TRUE(sss.empty());
        for (uint64_t k=0; k<1000000; ++k) {
            ASSERT_EQ(exp.size(), sss.size());
            uint64_t value = k+dice();
            if (exp.empty()) {
                exp.push(value);
                sss.push(value);
            } else {
                ASSERT_EQ(exp.top(), sss.top());
                if (exp.top() >= value) {
                    exp.pop();
                    sss.pop();
                } else {
                    exp.push(value);
                    sss.push(value);
                }
            }
        }
    }
}

//! Test Serialize and Load
TEST_F(sorted_stack_support_test, serialize_and_load)
{
    std::string file_name = temp_dir+"/sorted_stack_support";
    std::stack<uint64_t> exp;
    {
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, 10);
        auto dice = bind(distribution, rng);
        sdsl::sorted_stack_support sss(1000000+10);
        for (uint64_t k=0; k<1000000; ++k) {
            uint64_t value = k+dice();
            if (exp.empty() or exp.top() < value) {
                exp.push(value);
                sss.push(value);
            }
        }
        sdsl::store_to_file(sss, file_name);
    }
    {
        sdsl::sorted_stack_support sss(0);
        sdsl::load_from_file(sss, file_name);
        compare_stacks(exp, sss);
    }
    sdsl::remove(file_name);
}



}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        std::cout << "Usage: " << argv[0] << " tmp_dir" << std::endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    temp_dir = argv[1];
    return RUN_ALL_TESTS();
}
