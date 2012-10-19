#include "sdsl/suffixarrays.hpp"
#include "sdsl/config.hpp" // for CMAKE_SOURCE_DIR
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <algorithm> // for std::min

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
std::vector<sdsl::tMSS>  test_cases_file_map;

template<class T>
class CsaAsciiTest : public ::testing::Test
{
    protected:

        CsaAsciiTest() {
            // You can do set-up work for each test here.
        }

        virtual ~CsaAsciiTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            string test_cases_dir = string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/test_cases";
            tmp_dir = string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/tmp/";
            test_cases.push_back(test_cases_dir + "/crafted/100a.txt");
            test_cases.push_back(test_cases_dir + "/small/faust.txt");
            test_cases.push_back(test_cases_dir + "/small/zarathustra.txt");
            test_cases.push_back(test_cases_dir + "/crafted/empty.txt");
            tmp_file = "tmp_csa_ascii_test_" + sdsl::util::to_string(sdsl::util::get_pid()) + "_";
			if ( test_cases_file_map.size() == 0 ){
				test_cases_file_map.resize(test_cases.size());
			}
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        std::vector<std::string> test_cases;
        std::string tmp_file;
        std::string tmp_dir;

        template<class Csa>
        std::string get_tmp_file_name(const Csa& csa, size_type i) {
            std::locale loc;                 // the "C" locale
            const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc);
            std::string name = sdsl::util::demangle2(typeid(Csa).name());
            uint64_t myhash = coll.hash(name.data(),name.data()+name.length());
            return tmp_file + sdsl::util::to_string(myhash) + "_" + sdsl::util::basename(test_cases[i].c_str());
        }

        template<class Csa>
        bool load_csa(Csa& csa, size_type i) {
            return sdsl::util::load_from_file(csa, get_tmp_file_name(csa, i).c_str());
        }
};


using testing::Types;

typedef Types<  sdsl::csa_wt<>,
				sdsl::csa_sada<>,
				sdsl::csa_uncompressed	
		     > Implementations;

TYPED_TEST_CASE(CsaAsciiTest, Implementations);

TYPED_TEST(CsaAsciiTest, CreateAndStoreTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		TypeParam csa;
        construct_csa( this->test_cases[i].c_str(), csa, test_cases_file_map[i], false, 
				       this->tmp_dir, sdsl::util::basename(this->test_cases[i].c_str()) );
        bool success = sdsl::util::store_to_file(csa, this->get_tmp_file_name(csa, i).c_str());
        ASSERT_EQ(success, true);
    }
}

//! Test access methods
TYPED_TEST(CsaAsciiTest, Sigma)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(this->load_csa(csa, i), true);
        char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), text);
        ASSERT_EQ(n, csa.size());
        sdsl::bit_vector occur(256, 0);
        uint16_t sigma = 0;
        for (size_type j=0; j<n; ++j) {
            if (!occur[(unsigned char)text[j]]) {
                occur[(unsigned char)text[j]] = 1;
                ++sigma;
            }
        }
        ASSERT_EQ(sigma, csa.sigma);
        delete [] text;
    }
}

//! Test access methods
TYPED_TEST(CsaAsciiTest, Access) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(this->load_csa(csa, i), true);
		sdsl::int_vector<> sa;
		sdsl::util::load_from_file(sa, test_cases_file_map[i]["sa"].c_str());
        size_type n = sa.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(sa[j], csa[j])<<" j="<<j;
        }
    }
}

/*
template<class tWt>
void test_rank(const tWt& wt, const unsigned char* text, size_type n)
{
    std::vector<size_type> cnt(256, 0);
    ASSERT_EQ(n, wt.size());
    for (size_type j=0; j < wt.size(); ++j) {
        cnt[text[j]]++;
        ASSERT_EQ(cnt[text[j]], wt.rank(j+1, text[j]))<< " j = "<<j<<" text[j]"<<text[j];
    }
    // Do random queries for all characters that do not occur in the string
    for (size_type j=0; j<cnt.size(); ++j) {
        if (cnt[j] == 0) {
            for (size_type k=0; k<1000; ++k) {
                size_type pos = rand()%(wt.size()+1);
                ASSERT_EQ((size_type)0, wt.rank(pos, (unsigned char)j))<<" pos="<<pos;
            }
        }
    }
    // Test rank(size(), c) for each character c
    for (size_type c=0; c < 256; ++c) {
        ASSERT_EQ(cnt[c], wt.rank(wt.size(), (unsigned char)c))<<" c="<<c;
    }
}


//! Test rank methods
TYPED_TEST(CsaAsciiTest, Rank)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        unsigned char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), (char*&)text)-1;
        ::test_rank(wt, text, n);
        delete [] text;
    }
}


//! Test select methods
TYPED_TEST(CsaAsciiTest, Select)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        unsigned char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), (char*&)text)-1;
        std::vector<size_type> cnt(256, 0);
        ASSERT_EQ(n, wt.size());
        for (size_type j=0; j<n; ++j) {
            cnt[text[j]]++;
            ASSERT_EQ(j, wt.select(cnt[text[j]], text[j]))<< " j = "<<j<<" text[j]"<<text[j];
        }
        delete [] text;
    }
}
*/

//! Test access after swap
TYPED_TEST(CsaAsciiTest, SwapTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa1;
        ASSERT_EQ(this->load_csa(csa1, i), true);
        TypeParam csa2;
        csa1.swap(csa2);
		sdsl::int_vector<> sa;
		sdsl::util::load_from_file(sa, test_cases_file_map[i]["sa"].c_str());
        size_type n = sa.size();
        ASSERT_EQ(n, csa2.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(csa2[j], (typename TypeParam::value_type)sa[j]);
        }
    }
}


TYPED_TEST(CsaAsciiTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        std::remove(this->get_tmp_file_name(csa, i).c_str());
		sdsl::util::delete_all_files(test_cases_file_map[i]);	
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
