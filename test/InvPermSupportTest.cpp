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

// The fixture for testing class int_vector.
class InvPermSupportTest : public ::testing::Test
{
    protected:

        InvPermSupportTest() {}

        virtual ~InvPermSupportTest() {}

        virtual void SetUp() {
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

//! Test Constructors
TEST_F(InvPermSupportTest, Constructors)
{
    for (size_t i=0; i<perms.size(); ++i) {
        sdsl::inv_perm_support<> ips(&perms[i]);
        ASSERT_EQ(inv_perms[i].size(), ips.size());
        for (size_t j=0; j<ips.size(); ++j) {
            ASSERT_EQ(inv_perms[i][j], ips[j]);
        }
    }
}

TEST_F(InvPermSupportTest, Swap)
{
    for (size_type i=0; i < perms.size(); ++i) {
        sdsl::inv_perm_support<> ips(&perms[i]);
        {
            sdsl::inv_perm_support<> tmp;
            ASSERT_EQ((size_type)0, tmp.size());
            sdsl::util::swap_support(tmp, ips,
                                     &perms[i], (const sdsl::int_vector<>*)nullptr);
            ASSERT_EQ((size_type)0, ips.size());
            ASSERT_EQ(inv_perms[i].size(), tmp.size());
            for (size_type j=0; j < tmp.size(); ++j) {
                ASSERT_EQ(inv_perms[i][j], tmp[j]);
            }
        }
    }
}

TEST_F(InvPermSupportTest, SerializeAndLoad)
{
    for (size_type i=0; i < perms.size(); ++i) {
        std::string file_name = "tmp/inv_perm_support";
        {
            sdsl::inv_perm_support<> ips(&perms[i]);
            sdsl::store_to_file(ips, file_name);
        }
        sdsl::inv_perm_support<> ips2;
        sdsl::load_from_file(ips2, file_name);
        ips2.set_vector(&perms[i]);
        ASSERT_EQ(inv_perms[i].size(), ips2.size());
        for (size_t j=0; j<ips2.size(); ++j) {
            ASSERT_EQ(inv_perms[i][j], ips2[j]);
        }
        sdsl::remove(file_name);
    }
}

TEST_F(InvPermSupportTest, IteratorTest)
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
    return RUN_ALL_TESTS();
}
