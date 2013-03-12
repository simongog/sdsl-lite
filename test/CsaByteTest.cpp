#include "sdsl/suffix_arrays.hpp"
#include "sdsl/config.hpp" // for CMAKE_SOURCE_DIR
#include "sdsl/testutils.hpp"
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
class CsaByteTest : public ::testing::Test {
	protected:
	CsaByteTest() { }

	virtual ~CsaByteTest() { }

	virtual void SetUp() {
		tmp_dir = std::string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/tmp/";	
		std::string prefix		= std::string(SDSL_XSTR(CMAKE_SOURCE_DIR))+"/test";
		std::string config_file = prefix + "/CsaByteTest.config";
		std::string tc_prefix	= prefix + "/test_cases";
		test_cases = sdsl::paths_from_config_file(config_file, tc_prefix.c_str());
		tmp_file = "csa_ascii_test_" + util::to_string(util::pid()) + "_";
		if ( test_cases_file_map.size() == 0 ){
			test_cases_file_map.resize(test_cases.size());
		}
	}

	virtual void TearDown() {}

	std::vector<std::string> test_cases;
	std::string tmp_file;
	std::string tmp_dir;

	template<class Csa>
	std::string get_tmp_file_name(const Csa& csa, size_type i) {
		return tmp_dir+ tmp_file + util::class_to_hash(csa) + "_" + util::basename(test_cases[i]);
	}

	template<class Csa>
	bool load_csa(Csa& csa, size_type i) {
		return util::load_from_file(csa, get_tmp_file_name(csa, i));
	}
};


using testing::Types;

typedef Types<  csa_wt<>,
				csa_wt<wt_huff<>, 8, 16, text_order_sa_sampling<> >,
				csa_wt<wt_huff<>, 8, 16, sa_order_sa_sampling<> >,
				csa_wt<wt_huff<>, 8, 16, sa_order_sa_sampling<>, int_vector<>, 
				             succinct_byte_alphabet_strategy<bit_vector, rank_support_v<>, select_support_mcl<> > >,
				csa_wt<wt_huff<>, 8, 16, sa_order_sa_sampling<>, int_vector<>, 
				             succinct_byte_alphabet_strategy<> >,
				csa_sada<>,
				csa_bitcompressed<> 
		     > Implementations;

TYPED_TEST_CASE(CsaByteTest, Implementations);

TYPED_TEST(CsaByteTest, CreateAndStoreTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		TypeParam csa;
		cache_config config(false, this->tmp_dir, util::basename(this->test_cases[i]));
		construct(csa, this->test_cases[i], config, 1);
		test_cases_file_map[i] = config.file_map;
        bool success = util::store_to_file(csa, this->get_tmp_file_name(csa, i));
        ASSERT_EQ(true, success);
    }
}

//! Test sigma member
TYPED_TEST(CsaByteTest, Sigma) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(true, this->load_csa(csa, i));
		int_vector<8> text;
		ASSERT_EQ(true, util::load_vector_from_file(text, this->test_cases[i], 1));
		text.resize(text.size()+1);
		text[text.size()-1] = 0; // add 0-character to the end
        ASSERT_EQ(text.size(), csa.size());
        bit_vector occur(256, 0);
        uint16_t sigma = 0;
        for (size_type j=0; j<text.size(); ++j) {
            if (!occur[text[j]]) {
                occur[text[j]] = 1;
                ++sigma;
            }
        }
        ASSERT_EQ(sigma, csa.sigma);
    }
}

//! Test suffix array access methods
TYPED_TEST(CsaByteTest, SaAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(true, this->load_csa(csa, i));
		int_vector<> sa;
		util::load_from_file(sa, test_cases_file_map[i][constants::KEY_SA]);
        size_type n = sa.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(sa[j], csa[j])<<" j="<<j;
        }
    }
}

//! Test inverse suffix access methods
TYPED_TEST(CsaByteTest, IsaAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(true, this->load_csa(csa, i));
		int_vector<> isa;
		size_type n = 0;
		{
			int_vector<> sa;
			util::load_from_file(sa, test_cases_file_map[i][constants::KEY_SA]);
			n = sa.size();
			ASSERT_EQ(n, csa.size());
			isa = sa;
			for (size_type j=0; j<n; ++j){
				isa[sa[j]] = j;
			}
		}
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(isa[j], csa(j))<<" j="<<j;
        }
    }
}

//! Test text access methods
TYPED_TEST(CsaByteTest, TextAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		if ( test_cases_file_map[i].find(constants::KEY_TEXT) != test_cases_file_map[i].end() ){
			TypeParam csa;
			ASSERT_EQ(true, this->load_csa(csa, i));
			int_vector<8> text;
			util::load_from_file(text, test_cases_file_map[i][constants::KEY_TEXT]);
			size_type n = text.size();
			ASSERT_EQ(n, csa.size());
			for (size_type j=0; j<n; ++j) {
				ASSERT_EQ(text[j], csa.text[j])<<" j="<<j;
			}
		}
    }
}

//! Test Burrows-Wheeler access methods
TYPED_TEST(CsaByteTest, BwtAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		if ( test_cases_file_map[i].find(constants::KEY_BWT) != test_cases_file_map[i].end() ){
			TypeParam csa;
			ASSERT_EQ(true, this->load_csa(csa, i));
			int_vector<8> bwt;
			util::load_from_file(bwt, test_cases_file_map[i][constants::KEY_BWT]);
			size_type n = bwt.size();
			ASSERT_EQ(n, csa.size());
			for (size_type j=0; j<n; ++j) {
				ASSERT_EQ(bwt[j], csa.bwt[j])<<" j="<<j;
			}
		}
    }
}

//! Test Psi access methods
TYPED_TEST(CsaByteTest, PsiAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		if ( test_cases_file_map[i].find(constants::KEY_PSI) != test_cases_file_map[i].end() ){
			TypeParam csa;
			ASSERT_EQ(true, this->load_csa(csa, i));
			int_vector<> psi;
			util::load_from_file(psi, test_cases_file_map[i][constants::KEY_PSI]);
			size_type n = psi.size();
			ASSERT_EQ(n, csa.size());
			for (size_type j=0; j<n; ++j) {
				ASSERT_EQ(psi[j], csa.psi[j])<<" j="<<j;
			}
		}
    }
}

//! Test if Psi[LF[i]]=i 
TYPED_TEST(CsaByteTest, PsiLFAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		TypeParam csa;
		ASSERT_EQ(true, this->load_csa(csa, i));
		for (size_type j=0; j<csa.size(); ++j) {
			size_type lf = csa.psi(j);
			ASSERT_TRUE( lf >= 0 );
			ASSERT_TRUE( lf < csa.size() );
			ASSERT_EQ(j, csa.psi[lf])<<" j="<<j;
		}
    }
}


//! Test access after swap
TYPED_TEST(CsaByteTest, SwapTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa1;
        ASSERT_EQ(true, this->load_csa(csa1, i));
        TypeParam csa2;
        csa1.swap(csa2);
		int_vector<> sa;
		util::load_from_file(sa, test_cases_file_map[i][constants::KEY_SA]);
        size_type n = sa.size();
        ASSERT_EQ(n, csa2.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ((typename TypeParam::value_type)sa[j], csa2[j]);
        }
    }
}


TYPED_TEST(CsaByteTest, DeleteTest) {
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
