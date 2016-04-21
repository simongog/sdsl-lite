#include "sdsl/inv_perm_support.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <random>
#include <algorithm>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::int_vector<>::value_type value_type;

std::string temp_dir;

// The fixture for testing class int_vector.
class inv_perm_support_test : public ::testing::Test
{
    protected:

        inv_perm_support_test() {}

        virtual ~inv_perm_support_test() {}

        virtual void SetUp()
        {
            for (size_t z=0; z < (1ULL<<20); z=(z+1)*2) {
                sdsl::int_vector<> iv(z);
                sdsl::util::set_to_id(iv);
                random_shuffle(iv.begin(), iv.end());
                perms.emplace_back(iv);
            }
            for (size_t i=0; i<perms.size(); ++i) {
                auto iiv = perms[i];
                for (size_t j=0; j<perms[i].size(); ++j) {
                    iiv[perms[i][j]] = j;
                }
                inv_perms.emplace_back(iiv);
            }
        }

        virtual void TearDown() {}

        std::vector<sdsl::int_vector<>> perms;
        std::vector<sdsl::int_vector<>> inv_perms;
};

void compare(const sdsl::int_vector<>& inv_perm, const sdsl::inv_perm_support<>& ips)
{
    ASSERT_EQ(inv_perm.size(), ips.size());
    for (size_t j=0; j<ips.size(); ++j) {
        ASSERT_EQ(inv_perm[j], ips[j]);
    }
}

//! Test Constructors
TEST_F(inv_perm_support_test, constructors)
{
    static_assert(sdsl::util::is_regular<sdsl::inv_perm_support<>>::value, "Type is not regular");
    for (size_t i=0; i<perms.size(); ++i) {
        // Constructor
        sdsl::inv_perm_support<> ips(&perms[i]);
        compare(inv_perms[i], ips);

        // Copy-constructor
        sdsl::inv_perm_support<> ips2(ips);
        compare(inv_perms[i], ips2);

        // Move-constructor
        sdsl::inv_perm_support<> ips3(std::move(ips2));
        compare(inv_perms[i], ips3);

        // Copy-assign
        sdsl::inv_perm_support<> ips4;
        ips4 = ips;
        compare(inv_perms[i], ips4);

        // Move-assign
        sdsl::inv_perm_support<> ips5;
        ips5 = std::move(ips);
        compare(inv_perms[i], ips5);
    }
}

TEST_F(inv_perm_support_test, swap)
{
    for (size_type i=0; i < perms.size(); ++i) {
        sdsl::inv_perm_support<> ips(&perms[i]);
        {
            sdsl::inv_perm_support<> tmp;
            ASSERT_EQ((size_type)0, tmp.size());
            sdsl::util::swap_support(tmp, ips,
                                     &perms[i], (const sdsl::int_vector<>*)nullptr);
            ASSERT_EQ((size_type)0, ips.size());
            compare(inv_perms[i], tmp);
        }
    }
}

TEST_F(inv_perm_support_test, serialize_and_load)
{
    for (size_type i=0; i < perms.size(); ++i) {
        std::string file_name = temp_dir+"/inv_perm_support";
        {
            sdsl::inv_perm_support<> ips(&perms[i]);
            sdsl::store_to_file(ips, file_name);
        }
        sdsl::inv_perm_support<> ips2;
        sdsl::load_from_file(ips2, file_name);
        ips2.set_vector(&perms[i]);
        compare(inv_perms[i], ips2);
        sdsl::remove(file_name);
    }
}

TEST_F(inv_perm_support_test, iterator_test)
{
    for (size_type i=0; i < perms.size(); ++i) {
        sdsl::inv_perm_support<> ips(&perms[i]);
        auto iv_it  = inv_perms[i].begin();
        auto ips_it = ips.begin();
        while (ips_it != ips.end()) {
            ASSERT_TRUE(iv_it != inv_perms[i].end());
            ASSERT_EQ(*iv_it, *ips_it);
            ++iv_it;
            ++ips_it;
        }
    }
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
