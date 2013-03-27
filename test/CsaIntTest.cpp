#include "sdsl/suffix_arrays.hpp"
#include "sdsl/config.hpp" // for CMAKE_SOURCE_DIR
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>
#include <algorithm> // for std::min

namespace
{

using namespace sdsl;

typedef int_vector<>::size_type size_type;
std::vector<tMSS>  test_cases_file_map;

template<class T>
class CsaIntTest : public ::testing::Test
{
    protected:
        CsaIntTest() { }

        virtual ~CsaIntTest() { }

        virtual void SetUp() {
            tmp_dir = std::string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/tmp/";

            std::string prefix		= std::string(SDSL_XSTR(CMAKE_SOURCE_DIR))+"/test";
            std::string config_file = prefix + "/CsaIntTest.config";
            std::string tc_prefix	= prefix + "/test_cases";
            std::vector<std::string> read_test_cases;
            read_test_cases = sdsl::paths_from_config_file(config_file, tc_prefix.c_str());

            for (size_t i=0; i < read_test_cases.size(); ++i) {
                std::string tc = read_test_cases[i];
                size_t extentsion_len = 4;
                bool valid = true;
                if ((valid = (tc.size() >= extentsion_len))) {
                    std::string extension = tc.substr(tc.size()-extentsion_len, tc.size()-1);
                    if ((valid = (extension == ".int" or extension == ".txt"))) {
                        test_cases.push_back(tc);
                        if (extension == ".int") {
                            num_bytes.push_back(8);
                        } else if (extension == ".txt") {
                            num_bytes.push_back(1);
                        }
                    }
                }
                if (!valid) {
                    std::cerr << "Test case `"<<read_test_cases[i]<<"` was ignored." << std::endl;
                    std::cerr << "It did not have the extension `.int` or `.txt`" << std::endl;
                }
            }
            tmp_file = "tmp_csa_int_test_" + util::to_string(util::pid()) + "_";
            if (test_cases_file_map.size() == 0) {
                test_cases_file_map.resize(test_cases.size());
            }
        }

        virtual void TearDown() { }

        std::vector<std::string> test_cases;
        std::vector<uint8_t> num_bytes;
        std::string tmp_file;
        std::string tmp_dir;

        template<class Csa>
        std::string get_tmp_file_name(const Csa& csa, size_type i) {
            return tmp_file + util::class_to_hash(csa) + "_" + util::basename(test_cases[i]);
        }

        template<class Csa>
        bool load_csa(Csa& csa, size_type i) {
            return load_from_file(csa, get_tmp_file_name(csa, i));
        }
};


using testing::Types;

typedef Types<  csa_wt<wt_int<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> >,
        csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> >,
        csa_bitcompressed<int_alphabet_strategy<> >,
        csa_wt<wt_int<rrr_vector<63> >, 8, 8, sa_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> >,
        csa_wt<wt_int<>, 32, 32, text_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> >,
        csa_sada<enc_vector<>, 32, 32, text_order_sa_sampling<>, int_vector<>, int_alphabet_strategy<> >
        > Implementations;

TYPED_TEST_CASE(CsaIntTest, Implementations);

TYPED_TEST(CsaIntTest, CreateAndStoreTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        cache_config config(false, this->tmp_dir, util::basename(this->test_cases[i]));
        construct(csa, this->test_cases[i], config, this->num_bytes[i]);
        test_cases_file_map[i] = config.file_map;
        bool success = store_to_file(csa, this->get_tmp_file_name(csa, i));
        ASSERT_EQ(true, success);
    }
}

//! Test access methods
TYPED_TEST(CsaIntTest, Sigma)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(true, this->load_csa(csa, i));
        int_vector<> text;
        load_vector_from_file(text, this->test_cases[i], this->num_bytes[i]);
        text.resize(text.size()+1);
        text[text.size()-1] = 0;    // add 0-character at the end
        size_type n = text.size();
        ASSERT_EQ(n, csa.size());
        std::set<uint64_t> occur;
        size_type sigma = 0;
        for (size_type j=0; j<n; ++j) {
            if (occur.end() == occur.find(text[j])) {
                occur.insert(text[j]);
                ++sigma;
            }
        }
        ASSERT_EQ(sigma, csa.sigma);
    }
}

//! Test suffix array access methods
TYPED_TEST(CsaIntTest, SaAccess)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(true, this->load_csa(csa, i));
        int_vector<> sa;
        load_from_file(sa, test_cases_file_map[i][constants::KEY_SA]);
        size_type n = sa.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(sa[j], csa[j])<<" j="<<j;
        }
    }
}


//! Test inverse suffix access methods
TYPED_TEST(CsaIntTest, IsaAccess)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(true, this->load_csa(csa, i));
        int_vector<> isa;
        size_type n = 0;
        {
            int_vector<> sa;
            load_from_file(sa, test_cases_file_map[i][constants::KEY_SA]);
            n = sa.size();
            ASSERT_EQ(n, csa.size());
            isa = sa;
            for (size_type j=0; j<n; ++j) {
                isa[sa[j]] = j;    // calculate inverse suffix array
            }
        }
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(isa[j], csa(j))<<" j="<<j;
        }
    }
}

//! Test Burrows-Wheeler access methods
TYPED_TEST(CsaIntTest, BwtAccess)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        if (test_cases_file_map[i].end() != test_cases_file_map[i].find(constants::KEY_BWT_INT)) {
            TypeParam csa;
            ASSERT_EQ(true, this->load_csa(csa, i));
            int_vector<> bwt;
            load_from_file(bwt, test_cases_file_map[i][constants::KEY_BWT_INT]);
            size_type n = bwt.size();
            ASSERT_EQ(n, csa.size());
            for (size_type j=0; j<n; ++j) {
                ASSERT_EQ(bwt[j], csa.bwt[j])<<" j="<<j;
            }
        }
    }
}

//! Test text access methods
TYPED_TEST(CsaIntTest, TextAccess)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        if (test_cases_file_map[i].end() != test_cases_file_map[i].find(constants::KEY_TEXT_INT)) {
            TypeParam csa;
            ASSERT_EQ(true, this->load_csa(csa, i));
            int_vector<> text;
            load_from_file(text, test_cases_file_map[i][constants::KEY_TEXT_INT]);
            size_type n = text.size();
            ASSERT_EQ(n, csa.size());
            for (size_type j=0; j<n; ++j) {
                ASSERT_EQ(text[j], csa.text[j])<<" j="<<j;
            }
        }
    }
}

//! Test Psi access methods
TYPED_TEST(CsaIntTest, PsiAccess)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        if (test_cases_file_map[i].end() != test_cases_file_map[i].find(constants::KEY_PSI)) {
            TypeParam csa;
            ASSERT_EQ(true, this->load_csa(csa, i));
            int_vector<> psi;
            load_from_file(psi, test_cases_file_map[i][constants::KEY_PSI]);
            size_type n = psi.size();
            ASSERT_EQ(n, csa.size());
            for (size_type j=0; j<n; ++j) {
                ASSERT_EQ(psi[j], csa.psi[j])<<" j="<<j;
            }
        }
    }
}

//! Test if Psi[LF[i]]=i
TYPED_TEST(CsaIntTest, PsiLFAccess)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(true, this->load_csa(csa, i));
        for (size_type j=0; j<csa.size(); ++j) {
            size_type lf = csa.psi(j);
            ASSERT_TRUE(lf >= 0);
            ASSERT_TRUE(lf < csa.size());
            ASSERT_EQ(j, csa.psi[lf])<<" j="<<j;
        }
    }
}


//! Test access after swap
TYPED_TEST(CsaIntTest, SwapTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa1;
        ASSERT_EQ(true, this->load_csa(csa1, i));
        TypeParam csa2;
        csa1.swap(csa2);
        int_vector<> sa;
        load_from_file(sa, test_cases_file_map[i][constants::KEY_SA]);
        size_type n = sa.size();
        ASSERT_EQ(n, csa2.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ((typename TypeParam::value_type)sa[j], csa2[j]);
        }
    }
}

TYPED_TEST(CsaIntTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        std::remove(this->get_tmp_file_name(csa, i).c_str());
        util::delete_all_files(test_cases_file_map[i]);
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
