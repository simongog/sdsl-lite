#include "sdsl/wt_int.hpp"
#include "sdsl/construct.hpp"
#include "sdsl/util.hpp"
#include "sdsl/rrr_vector.hpp"
#include "sdsl/sd_vector.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <map>

namespace
{

using namespace sdsl;

typedef int_vector<>::size_type size_type;
typedef std::map<int_vector<>::value_type,size_type> tMII;

template<class T>
class WtIntTest : public ::testing::Test
{
    protected:

        WtIntTest() { }

        virtual ~WtIntTest() { }

        virtual void SetUp() {
            std::string tmp_dir = std::string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/tmp/";
            tmp_file = tmp_dir+"wt_int_test_" + util::to_string(util::pid()) + "_";
            std::string prefix		= std::string(SDSL_XSTR(CMAKE_SOURCE_DIR))+"/test";
            std::string config_file = prefix + "/WtIntTest.config";
            std::string tc_prefix	= prefix + "/test_cases";
            std::vector<std::string> read_test_cases;
            test_cases = sdsl::paths_from_config_file(config_file, tc_prefix.c_str());
        }

        virtual void TearDown() { }

        template<class Wt>
        std::string get_tmp_file_name(const Wt& wt, size_type i) {
            return tmp_file + util::class_to_hash(wt) + "_" + util::to_string(i);
        }

        template<class Wt>
        bool load_wt(Wt& wt, size_type i) {
            return load_from_file(wt, get_tmp_file_name(wt, i));
        }
        std::string tmp_file;
        std::vector<std::string> test_cases;
};

using testing::Types;

typedef Types<
wt_int<rrr_vector<15> >,
       wt_int<>,
       wt_int<rrr_vector<63> >
       > Implementations;

TYPED_TEST_CASE(WtIntTest, Implementations);

//! Test the parametrized constructor
TYPED_TEST(WtIntTest, Constructor)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        int_vector<> iv;
        load_from_file(iv, this->test_cases[i]);
        double iv_size = size_in_mega_bytes(iv);
        std::cout << "tc = " << this->test_cases[i] << std::endl;
        {
            TypeParam wt;
            sdsl::construct(wt, this->test_cases[i]);
            std::cout << "compression = " << size_in_mega_bytes(wt)/iv_size << std::endl;
            ASSERT_EQ(iv.size(), wt.size());
            for (size_type j=0; j < iv.size(); ++j) {
                ASSERT_EQ(iv[j], wt[j])<<j;
            }
            ASSERT_TRUE(store_to_file(wt, this->get_tmp_file_name(wt, i)));
        }
        {
            int_vector_file_buffer<> iv_buf(this->test_cases[i]);
            TypeParam wt(iv_buf, 0);
            ASSERT_EQ((size_type)0,  wt.size());
        }
        {
            int_vector_file_buffer<> iv_buf(this->test_cases[i]);
            size_type len = (iv.size() >= 6) ? 6 : iv.size();
            TypeParam wt(iv_buf, len);
            ASSERT_EQ(len, wt.size());
            for (size_type j=0; j < len; ++j) {
                ASSERT_EQ(iv[j], wt[j])<<j;
            }
        }
    }
}

//! Test loading and accessing the wavelet tree
TYPED_TEST(WtIntTest, LoadAndAccess)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        int_vector<> iv;
        load_from_file(iv, this->test_cases[i]);
        std::string tmp_file_name = this->tmp_file+util::to_string(i);
        TypeParam wt;
        ASSERT_TRUE(load_from_file(wt, this->get_tmp_file_name(wt, i)));
        ASSERT_EQ(iv.size(), wt.size());
        for (size_type j=0; j < iv.size(); ++j) {
            ASSERT_EQ(iv[j], wt[j])<<j;
        }
    }
}

//! Test the load method and rank method
TYPED_TEST(WtIntTest, LoadAndRank)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        int_vector<> iv;
        load_from_file(iv, this->test_cases[i]);
        std::string tmp_file_name = this->tmp_file+util::to_string(i);
        TypeParam wt;
        ASSERT_TRUE(load_from_file(wt, this->get_tmp_file_name(wt, i)));
        ASSERT_EQ(iv.size(), wt.size());
        tMII check_rank;
        for (size_type j=0; j < iv.size(); ++j) {
            ASSERT_EQ(wt.rank(j, iv[j]), check_rank[iv[j]]);
            check_rank[iv[j]]++;
        }
        for (tMII::const_iterator it=check_rank.begin(); it!=check_rank.end(); ++it) {
            ASSERT_EQ(wt.rank(wt.size(), it->first), it->second);
        }
    }
}

//! Test the load method and select method
TYPED_TEST(WtIntTest, LoadAndSelect)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        int_vector<> iv;
        load_from_file(iv, this->test_cases[i]);
        std::string tmp_file_name = this->tmp_file+util::to_string(i);
        TypeParam wt;
        ASSERT_TRUE(load_from_file(wt, this->get_tmp_file_name(wt, i)));
        ASSERT_EQ(iv.size(), wt.size());
        tMII count;
        for (size_type j=0; j < iv.size(); ++j) {
            count[iv[j]]++;
            ASSERT_EQ(j, wt.select(count[iv[j]], iv[j])) << "iv[j]=" << iv[j] << " "
                    << " j="<<j;
        }
    }
}

//! Test the load method and inverse_select method
TYPED_TEST(WtIntTest, LoadAndInverseSelect)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        int_vector<> iv;
        load_from_file(iv, this->test_cases[i]);
        std::string tmp_file_name = this->tmp_file+util::to_string(i);
        TypeParam wt;
        ASSERT_TRUE(load_from_file(wt, this->get_tmp_file_name(wt, i)));
        ASSERT_EQ(iv.size(), wt.size());
        tMII check_rank;
        for (size_type j=0; j < iv.size(); ++j) {
            typename TypeParam::value_type x;
            ASSERT_EQ(check_rank[iv[j]], wt.inverse_select(j, x));
            ASSERT_EQ(iv[j], x);
            check_rank[iv[j]]++;
        }
    }
}

TYPED_TEST(WtIntTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        std::string tmp_file_name = this->tmp_file+util::to_string(i);
        std::remove(tmp_file_name.c_str());
        TypeParam wt;
        std::remove(this->get_tmp_file_name(wt, i).c_str());
    }
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
