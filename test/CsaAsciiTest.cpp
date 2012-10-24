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
				sdsl::csa_wt<sdsl::wt_huff<>, 8, 16, sdsl::text_order_sa_sampling<> >,
				sdsl::csa_wt<sdsl::wt_huff<>, 8, 16, sdsl::sa_order_sa_sampling<> >,
				sdsl::csa_wt<sdsl::wt_huff<>, 8, 16, sdsl::sa_order_sa_sampling<>, sdsl::int_vector<>, 
				             sdsl::succinct_byte_alphabet_strategy<sdsl::bit_vector, sdsl::rank_support_v<>, sdsl::select_support_mcl<> > >,
				sdsl::csa_wt<sdsl::wt_huff<>, 8, 16, sdsl::sa_order_sa_sampling<>, sdsl::int_vector<>, 
				             sdsl::succinct_byte_alphabet_strategy<> >,
				sdsl::csa_sada<>,
				sdsl::csa_bitcompressed<> 
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

//! Test suffix array access methods
TYPED_TEST(CsaAsciiTest, SaAccess) {
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

//! Test inverse suffix access methods
TYPED_TEST(CsaAsciiTest, IsaAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam csa;
        ASSERT_EQ(this->load_csa(csa, i), true);
		sdsl::int_vector<> isa;
		size_type n = 0;
		{
			sdsl::int_vector<> sa;
			sdsl::util::load_from_file(sa, test_cases_file_map[i]["sa"].c_str());
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

//! Test Burrows-Wheeler access methods
TYPED_TEST(CsaAsciiTest, BwtAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		if ( test_cases_file_map[i].find("bwt") != test_cases_file_map[i].end() ){
			TypeParam csa;
			ASSERT_EQ(this->load_csa(csa, i), true);
			sdsl::int_vector<8> bwt;
			sdsl::util::load_from_file(bwt, test_cases_file_map[i]["bwt"].c_str());
			size_type n = bwt.size();
			ASSERT_EQ(n, csa.size());
			for (size_type j=0; j<n; ++j) {
				ASSERT_EQ(bwt[j], csa.bwt[j])<<" j="<<j;
			}
		}
    }
}

//! Test Psi access methods
TYPED_TEST(CsaAsciiTest, PsiAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		if ( test_cases_file_map[i].find("psi") != test_cases_file_map[i].end() ){
			TypeParam csa;
			ASSERT_EQ(this->load_csa(csa, i), true);
			sdsl::int_vector<> psi;
			sdsl::util::load_from_file(psi, test_cases_file_map[i]["psi"].c_str());
			size_type n = psi.size();
			ASSERT_EQ(n, csa.size());
			for (size_type j=0; j<n; ++j) {
				ASSERT_EQ(psi[j], csa.psi[j])<<" j="<<j;
			}
		}
    }
}

//! Test if Psi[LF[i]]=i 
TYPED_TEST(CsaAsciiTest, PsiLFAccess) {
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
