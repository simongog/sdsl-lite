#include "sdsl/wt_int.hpp"
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

typedef sdsl::int_vector<>::size_type size_type;
typedef std::map<sdsl::int_vector<>::value_type,size_type> tMII;

template<class T>
class WtIntTest : public ::testing::Test
{
    protected:

        WtIntTest() {
            // You can do set-up work for each test here.
        }

        virtual ~WtIntTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            tmp_file = "tmp_wt_int_test_" + sdsl::util::to_string(sdsl::util::get_pid()) + "_";
			test_cases.push_back( sdsl::int_vector<>() );
			test_cases.push_back( sdsl::int_vector<>(1023,0,1) );
			test_cases.push_back( sdsl::int_vector<>(100023,0,1) );
			test_cases.push_back( sdsl::int_vector<>(64,0,2) );
			sdsl::int_vector<> iv(1000000,0,18);
    		sdsl::util::set_random_bits(iv, 17);
			test_cases.push_back( iv );
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        template<class Wt>
        std::string get_tmp_file_name(const Wt& wt, size_type i) {
            return tmp_file + sdsl::util::class_to_hash(wt) + "_" + sdsl::util::to_string(i);
        }

        template<class Wt>
        bool load_wt(Wt& wt, size_type i) {
            return sdsl::util::load_from_file(wt, get_tmp_file_name(wt, i).c_str());
        }
        // Objects declared here can be used by all tests in the test case for Foo.
        std::string tmp_file;
        std::vector<sdsl::int_vector<> > test_cases;
};

using testing::Types;

typedef Types<
     		  sdsl::wt_int<sdsl::rrr_vector<15> >,
			  sdsl::wt_int<>,
     		  sdsl::wt_int<sdsl::rrr_vector<63> >
     		 > Implementations;

TYPED_TEST_CASE(WtIntTest, Implementations);

//! Test the parametrized constructor
TYPED_TEST(WtIntTest, Constructor) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		sdsl::int_vector<>& iv = this->test_cases[i];
		double iv_size = sdsl::util::get_size_in_mega_bytes(iv);
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
		ASSERT_TRUE( sdsl::util::store_to_file(iv, tmp_file_name.c_str()) );
		{
		sdsl::int_vector_file_buffer<> buf(tmp_file_name.c_str());
			TypeParam wt(buf, buf.int_vector_size);
			std::cout << "compression = " << sdsl::util::get_size_in_mega_bytes(wt)/iv_size << std::endl;
			ASSERT_EQ(iv.size(), wt.size());
			for (size_type j=0; j < iv.size(); ++j) {
				ASSERT_EQ(iv[j], wt[j])<<j;
			}
			ASSERT_TRUE( sdsl::util::store_to_file(wt, get_tmp_file_name(wt, i).c_str()) );
		}
		{
			sdsl::int_vector_file_buffer<> buf(tmp_file_name.c_str());
			TypeParam wt(buf, 0);
			ASSERT_EQ( wt.size(), (size_type)0 );
		}
		{
			sdsl::int_vector_file_buffer<> buf(tmp_file_name.c_str());
			size_type len = (iv.size() >= 6) ? 6 : iv.size(); 
			TypeParam wt(buf, len);
			ASSERT_EQ( wt.size(), len );
			for (size_type j=0; j < len; ++j) {
				ASSERT_EQ(iv[j], wt[j])<<j;
			}
		}
	}
}

//! Test loading and accessing the wavelet tree
TYPED_TEST(WtIntTest, LoadAndAccess) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		sdsl::int_vector<>& iv = this->test_cases[i];
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
		TypeParam wt;
		ASSERT_TRUE(sdsl::util::load_from_file(wt, get_tmp_file_name(wt, i).c_str()));
		ASSERT_EQ(iv.size(), wt.size());
		for (size_type j=0; j < iv.size(); ++j) {
			ASSERT_EQ(iv[j], wt[j])<<j;
		}
	}
}

//! Test the load method and rank method
TYPED_TEST(WtIntTest, LoadAndRank) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		sdsl::int_vector<>& iv = this->test_cases[i];
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
		TypeParam wt;
		ASSERT_TRUE(sdsl::util::load_from_file(wt, get_tmp_file_name(wt, i).c_str()));
		ASSERT_EQ(iv.size(), wt.size());
		tMII check_rank;
		for (size_type j=0; j < iv.size(); ++j) {
			ASSERT_EQ(wt.rank(j, iv[j]), check_rank[iv[j]]);
			check_rank[iv[j]]++;
		}
		for (tMII::const_iterator it=check_rank.begin(); it!=check_rank.end(); ++it){
			ASSERT_EQ(wt.rank(wt.size(), it->first), it->second);
		}
	}
}


//! Test the load method and select method
TYPED_TEST(WtIntTest, LoadAndSelect) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		sdsl::int_vector<>& iv = this->test_cases[i];
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
		TypeParam wt;
		ASSERT_TRUE(sdsl::util::load_from_file(wt, get_tmp_file_name(wt, i).c_str()));
		ASSERT_EQ(iv.size(), wt.size());
		tMII count;
		for (size_type j=0; j < iv.size(); ++j) {
			count[iv[j]]++;
			ASSERT_EQ(wt.select(count[iv[j]], iv[j]), j) << "iv[j]=" << iv[j] << " "
			         << " j="<<j; 
		}
	}
}


//! Test the load method and inverse_select method
TYPED_TEST(WtIntTest, LoadAndInverseSelect) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		sdsl::int_vector<>& iv = this->test_cases[i];
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
		TypeParam wt;
		ASSERT_TRUE(sdsl::util::load_from_file(wt, get_tmp_file_name(wt, i).c_str()));
		ASSERT_EQ(iv.size(), wt.size());
		tMII check_rank;
		for (size_type j=0; j < iv.size(); ++j) {
			typename TypeParam::value_type x;
			ASSERT_EQ(wt.inverse_select(j, x), check_rank[iv[j]]);
			ASSERT_EQ(x, iv[j]);
			check_rank[iv[j]]++;
		}
	}
}

TYPED_TEST(WtIntTest, DeleteTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
		std::string tmp_file_name = this->tmp_file+sdsl::util::to_string(i);
        std::remove(tmp_file_name.c_str());
		TypeParam wt;
		std::remove(get_tmp_file_name(wt, i).c_str());
    }
}


}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
