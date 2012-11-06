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
class CsaIntTest : public ::testing::Test
{
    protected:

        CsaIntTest() {
            // You can do set-up work for each test here.
        }

        virtual ~CsaIntTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            string test_cases_dir = string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/test_cases";
            tmp_dir = string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/tmp/";
            test_cases.push_back(test_cases_dir + "/small/keeper.int"); num_bytes.push_back(8);
            test_cases.push_back(test_cases_dir + "/small/moby.int"); num_bytes.push_back(8);
			test_cases.push_back(test_cases_dir + "/crafted/empty.txt"); num_bytes.push_back(8);
            test_cases.push_back(test_cases_dir + "/small/faust.txt"); num_bytes.push_back(1);
            test_cases.push_back(test_cases_dir + "/small/zarathustra.txt"); num_bytes.push_back(1);
            tmp_file = "tmp_csa_int_test_" + sdsl::util::to_string(sdsl::util::get_pid()) + "_";
			if ( test_cases_file_map.size() == 0 ){
				test_cases_file_map.resize(test_cases.size());
			}
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        std::vector<std::string> test_cases;
        std::vector<uint8_t> num_bytes;
        std::string tmp_file;
        std::string tmp_dir;

        template<class Csa>
        std::string get_tmp_file_name(const Csa& csa, size_type i) {
            return tmp_file + sdsl::util::class_to_hash(csa) + "_" + sdsl::util::basename(test_cases[i].c_str());
        }

        template<class Csa>
        bool load_csa(Csa& csa, size_type i) {
            return sdsl::util::load_from_file(csa, get_tmp_file_name(csa, i).c_str());
        }
};


using testing::Types;

typedef Types<  sdsl::csa_wt<sdsl::wt_int<>, 32, 32, sdsl::sa_order_sa_sampling<>, sdsl::int_vector<>, sdsl::int_alphabet_strategy<> >,
				sdsl::csa_sada<sdsl::enc_vector<>, 32, 32, sdsl::sa_order_sa_sampling<>, sdsl::int_vector<>, sdsl::int_alphabet_strategy<> >,
				sdsl::csa_bitcompressed<sdsl::int_alphabet_strategy<> >,
				sdsl::csa_wt<sdsl::wt_int<sdsl::rrr_vector<63> >, 8, 8, sdsl::sa_order_sa_sampling<>, sdsl::int_vector<>, sdsl::int_alphabet_strategy<> >,
				sdsl::csa_wt<sdsl::wt_int<>, 32, 32, sdsl::text_order_sa_sampling<>, sdsl::int_vector<>, sdsl::int_alphabet_strategy<> >,
				sdsl::csa_sada<sdsl::enc_vector<>, 32, 32, sdsl::text_order_sa_sampling<>, sdsl::int_vector<>, sdsl::int_alphabet_strategy<> >
		     > Implementations;

TYPED_TEST_CASE(CsaIntTest, Implementations);

TYPED_TEST(CsaIntTest, CreateAndStoreTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		TypeParam csa;
		// TODO: construct method has to me modified for this
		sdsl::cache_config config(false, this->tmp_dir, sdsl::util::basename(this->test_cases[i].c_str()));
		sdsl::construct(csa, this->test_cases[i].c_str(), config, this->num_bytes[i]);
		test_cases_file_map[i] = config.file_map;
        bool success = sdsl::util::store_to_file(csa, this->get_tmp_file_name(csa, i).c_str());
        ASSERT_EQ(success, true);
    }
}

//! Test access methods
TYPED_TEST(CsaIntTest, Sigma) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(this->load_csa(csa, i), true);
		sdsl::int_vector<> text;
        sdsl::util::load_vector_from_file(text, this->test_cases[i].c_str(), this->num_bytes[i]);
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
TYPED_TEST(CsaIntTest, SaAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(this->load_csa(csa, i), true);
		sdsl::int_vector<> sa;
		sdsl::util::load_from_file(sa, test_cases_file_map[i][sdsl::constants::KEY_SA].c_str());
        size_type n = sa.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(sa[j], csa[j])<<" j="<<j;
        }
    }
}


//! Test inverse suffix access methods
TYPED_TEST(CsaIntTest, IsaAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(this->load_csa(csa, i), true);
		sdsl::int_vector<> isa;
		size_type n = 0;
		{
			sdsl::int_vector<> sa;
			sdsl::util::load_from_file(sa, test_cases_file_map[i][sdsl::constants::KEY_SA].c_str());
			n = sa.size();
			ASSERT_EQ(n, csa.size());
			isa = sa;
			for (size_type j=0; j<n; ++j){ isa[sa[j]] = j; } // calculate inverse suffix array
		}
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(isa[j], csa(j))<<" j="<<j;
        }
    }
}

//! Test Burrows-Wheeler access methods
TYPED_TEST(CsaIntTest, BwtAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		if ( test_cases_file_map[i].end() != test_cases_file_map[i].find(sdsl::constants::KEY_BWT_INT)  ){
			TypeParam csa;
			ASSERT_EQ(this->load_csa(csa, i), true);
			sdsl::int_vector<> bwt;
			sdsl::util::load_from_file(bwt, test_cases_file_map[i][sdsl::constants::KEY_BWT_INT].c_str());
			size_type n = bwt.size();
			ASSERT_EQ(n, csa.size());
			for (size_type j=0; j<n; ++j) {
				ASSERT_EQ(bwt[j], csa.bwt[j])<<" j="<<j;
			}
		}
    }
}

//! Test Psi access methods
TYPED_TEST(CsaIntTest, PsiAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		if ( test_cases_file_map[i].end() != test_cases_file_map[i].find(sdsl::constants::KEY_PSI) ){
			TypeParam csa;
			ASSERT_EQ(this->load_csa(csa, i), true);
			sdsl::int_vector<> psi;
			sdsl::util::load_from_file(psi, test_cases_file_map[i][sdsl::constants::KEY_PSI].c_str());
			size_type n = psi.size();
			ASSERT_EQ(n, csa.size());
			for (size_type j=0; j<n; ++j) {
				ASSERT_EQ(psi[j], csa.psi[j])<<" j="<<j;
			}
		}
    }
}

//! Test if Psi[LF[i]]=i 
TYPED_TEST(CsaIntTest, PsiLFAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		TypeParam csa;
		ASSERT_EQ(this->load_csa(csa, i), true);
		for (size_type j=0; j<csa.size(); ++j) {
			size_type lf = csa.psi(j);
			ASSERT_TRUE( lf >= 0 );
			ASSERT_TRUE( lf < csa.size() );
			ASSERT_EQ(j, csa.psi[lf])<<" j="<<j;
		}
    }
}


//! Test access after swap
TYPED_TEST(CsaIntTest, SwapTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa1;
        ASSERT_EQ(this->load_csa(csa1, i), true);
        TypeParam csa2;
        csa1.swap(csa2);
		sdsl::int_vector<> sa;
		sdsl::util::load_from_file(sa, test_cases_file_map[i][sdsl::constants::KEY_SA].c_str());
        size_type n = sa.size();
        ASSERT_EQ(n, csa2.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(csa2[j], (typename TypeParam::value_type)sa[j]);
        }
    }
}

TYPED_TEST(CsaIntTest, DeleteTest) {
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
