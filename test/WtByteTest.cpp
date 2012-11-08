#include "sdsl/wavelet_trees.hpp"
#include "sdsl/rrr_vector.hpp"
#include "sdsl/bit_vector_interleaved.hpp"
#include "sdsl/util.hpp"
#include "sdsl/bwt_construct.hpp"
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

template<class T>
class WtByteTest : public ::testing::Test
{
    protected:

        WtByteTest() { }

        virtual ~WtByteTest() { }

        virtual void SetUp() {
			std::string test_cases_dir = std::string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/test_cases";
			std::string tmp_dir = std::string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/tmp/";
            test_cases.push_back(test_cases_dir + "/crafted/100a.txt");
            test_cases.push_back(test_cases_dir + "/small/faust.txt");
            test_cases.push_back(test_cases_dir + "/small/zarathustra.txt");
            test_cases.push_back(test_cases_dir + "/crafted/empty.txt");
            tmp_file = tmp_dir + "wt_ascii_test_" + util::to_string(util::get_pid()) + "_";
        }

        virtual void TearDown() { }

        std::vector<std::string> test_cases;
        std::string tmp_file;

        template<class Wt>
        std::string get_tmp_file_name(const Wt& wt, size_type i) {
            return tmp_file + util::class_to_hash(wt) + "_" + util::basename(test_cases[i].c_str());
        }

        template<class Wt>
        bool load_wt(Wt& wt, size_type i) {
            return util::load_from_file(wt, get_tmp_file_name(wt, i).c_str());
        }
};

using testing::Types;

typedef Types<
	 wt<unsigned char*, rrr_vector<63> >,
     wt<unsigned char*, bit_vector_interleaved<>  >,
     wt<unsigned char*, bit_vector>,
     wt_huff<bit_vector_interleaved<> >,
     wt_huff<bit_vector, rank_support_v<> >,
     wt_huff<bit_vector, rank_support_v5<> >,
     wt_huff<rrr_vector<63> >,
     wt_rlmn<>,
     wt_rlmn<bit_vector>,
     wt_rlg<>,
     wt_rlg8<>
     > Implementations;

TYPED_TEST_CASE(WtByteTest, Implementations);

TYPED_TEST(WtByteTest, CreateAndStoreTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        int_vector_file_buffer<8> text_buf;
        text_buf.load_from_plain(this->test_cases[i].c_str());
        size_type n = text_buf.int_vector_size;
        TypeParam wt(text_buf, n);
        bool success = util::store_to_file(wt, this->get_tmp_file_name(wt, i).c_str());
        ASSERT_EQ(success, true);
    }
}

//! Test access methods
TYPED_TEST(WtByteTest, Sigma) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        char* text = NULL;
        size_type n = file::read_text((this->test_cases[i]).c_str(), text)-1;
        ASSERT_EQ(n, wt.size());
        bit_vector occur(256, 0);
        uint16_t sigma = 0;
        for (size_type j=0; j<n; ++j) {
            if (!occur[(unsigned char)text[j]]) {
                occur[(unsigned char)text[j]] = 1;
                ++sigma;
            }
        }
        ASSERT_EQ(sigma, wt.sigma);
        delete [] text;
    }
}

//! Test access methods
TYPED_TEST(WtByteTest, Access) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        char* text = NULL;
        size_type n = file::read_text((this->test_cases[i]).c_str(), text)-1;
        ASSERT_EQ(n, wt.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ((typename TypeParam::value_type)text[j], wt[j])<<" j="<<j;
        }
        delete [] text;
    }
}

template<class tWt>
void test_rank(const tWt& wt, const unsigned char* text, size_type n) {
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
TYPED_TEST(WtByteTest, Rank) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        unsigned char* text = NULL;
        size_type n = file::read_text((this->test_cases[i]).c_str(), (char*&)text)-1;
        ::test_rank(wt, text, n);
        delete [] text;
    }
}

//! Test select methods
TYPED_TEST(WtByteTest, Select) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        unsigned char* text = NULL;
        size_type n = file::read_text((this->test_cases[i]).c_str(), (char*&)text)-1;
        std::vector<size_type> cnt(256, 0);
        ASSERT_EQ(n, wt.size());
        for (size_type j=0; j<n; ++j) {
            cnt[text[j]]++;
            ASSERT_EQ(j, wt.select(cnt[text[j]], text[j]))<< " j = "<<j<<" text[j]"<<text[j];
        }
        delete [] text;
    }
}


//! Test access after swap
TYPED_TEST(WtByteTest, SwapTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt1;
        ASSERT_EQ(this->load_wt(wt1, i), true);
        TypeParam wt2;
        wt1.swap(wt2);
        unsigned char* text = NULL;
        size_type n = file::read_text((this->test_cases[i]).c_str(), (char*&)text)-1;
        ASSERT_EQ(n, wt2.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(wt2[j], (typename TypeParam::value_type)text[j]);
        }
        delete [] text;
    }
}

TYPED_TEST(WtByteTest, CreatePartiallyTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        int_vector_file_buffer<8> text_buf;
        text_buf.load_from_plain(this->test_cases[i].c_str());
        unsigned char* text = NULL;
        size_type n = file::read_text((this->test_cases[i]).c_str(), (char*&)text)-1;
        n = std::min(n, (size_type)50);
        TypeParam wt(text_buf, n);
        ::test_rank(wt, text, n);
        delete [] text;
    }
}

TYPED_TEST(WtByteTest, DeleteTest) {
    for (size_t i=0; i< this->test_cases.size(); ++i) {
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
