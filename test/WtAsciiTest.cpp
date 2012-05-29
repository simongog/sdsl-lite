#include "sdsl/wt_huff.hpp"
#include "sdsl/wt_rlmn.hpp"
#include "sdsl/rrr_vector.hpp"
#include "sdsl/bit_vector_interleaved.hpp"
#include "sdsl/util.hpp"
#include "sdsl/csa_construct.hpp"
#include "sdsl/bwt_construct.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;

template<class T>
class WtAsciiTest : public ::testing::Test
{
    protected:

        WtAsciiTest() {
            // You can do set-up work for each test here.
        }

        virtual ~WtAsciiTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
            test_cases.push_back("test_cases/crafted/100a.txt");
            test_cases.push_back("test_cases/small/faust.txt");
//           test_cases.push_back("test_cases/small/zarathustra.txt");
            tmp_file = "tmp_wt_ascii_test_" + sdsl::util::to_string(sdsl::util::get_pid()) + "_";
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        std::vector<std::string> test_cases;
        std::string tmp_file;

        template<class Wt>
        std::string get_tmp_file_name(const Wt& wt, size_type i) {
            std::locale loc;                 // the "C" locale
            const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc);
            std::string name = sdsl::util::demangle2(typeid(Wt).name());
            uint64_t myhash = coll.hash(name.data(),name.data()+name.length());
            return tmp_file + sdsl::util::to_string(myhash) + "_" + sdsl::util::basename(test_cases[i].c_str());
        }

        template<class Wt>
        bool load_wt(Wt& wt, size_type i) {
            return sdsl::util::load_from_file(wt, get_tmp_file_name(wt, i).c_str());
        }
};

using testing::Types;

typedef Types<
sdsl::wt_huff<sdsl::bit_vector, sdsl::rank_support_v<> >,
     sdsl::wt_huff<sdsl::bit_vector, sdsl::rank_support_v5<> >,
     sdsl::wt_huff<sdsl::bit_vector_interleaved<> >,
     sdsl::wt_huff<sdsl::rrr_vector<63> >,
     sdsl::wt_rlmn<>,
     sdsl::wt_rlmn<sdsl::bit_vector>,
     sdsl::wt<unsigned char*, sdsl::bit_vector>,
     sdsl::wt<unsigned char*, sdsl::bit_vector_interleaved<>  >,
     sdsl::wt<unsigned char*, sdsl::rrr_vector<63> >
     > Implementations;

TYPED_TEST_CASE(WtAsciiTest, Implementations);

TYPED_TEST(WtAsciiTest, CreateAndStoreTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam dummy_wt;
        std::string tmp_text = "text_"+this->get_tmp_file_name(dummy_wt, i);
        char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), text);
        sdsl::util::store_to_file(sdsl::char_array_serialize_wrapper<>((unsigned char*)text, n), tmp_text.c_str());
        delete [] text;
        sdsl::int_vector_file_buffer<8> text_buf(tmp_text.c_str());
        TypeParam wt(text_buf, n);
        bool success = sdsl::util::store_to_file(wt, this->get_tmp_file_name(wt, i).c_str());
        std::remove(tmp_text.c_str());
        ASSERT_EQ(success, true);
    }
}

//! Test access methods
TYPED_TEST(WtAsciiTest, Access)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), text);
        ASSERT_EQ(wt.size() ,n);
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(wt[j], (typename TypeParam::value_type)text[j]);
        }
        delete [] text;
    }
}


//! Test rank methods
TYPED_TEST(WtAsciiTest, Rank)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        unsigned char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), (char*&)text);
        std::vector<size_type> cnt(256, 0);
        ASSERT_EQ(wt.size() ,n);
        for (size_type j=0; j<n; ++j) {
            cnt[text[j]]++;
            ASSERT_EQ(wt.rank(j+1, text[j]), cnt[text[j]])<< " j = "<<j<<" text[j]"<<text[j];
        }
        delete [] text;
    }
}


//! Test select methods
TYPED_TEST(WtAsciiTest, Select)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt;
        ASSERT_EQ(this->load_wt(wt, i), true);
        unsigned char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), (char*&)text);
        std::vector<size_type> cnt(256, 0);
        ASSERT_EQ(wt.size() ,n);
        for (size_type j=0; j<n; ++j) {
            cnt[text[j]]++;
            ASSERT_EQ(wt.select(cnt[text[j]], text[j]), j)<< " j = "<<j<<" text[j]"<<text[j];
        }
        delete [] text;
    }
}


//! Test access after swap
TYPED_TEST(WtAsciiTest, SwapTest)
{
    for (size_t i=0; i< this->test_cases.size(); ++i) {
        TypeParam wt1;
        ASSERT_EQ(this->load_wt(wt1, i), true);
        TypeParam wt2;
        wt1.swap(wt2);
        char* text = NULL;
        size_type n = sdsl::file::read_text((this->test_cases[i]).c_str(), text);
        ASSERT_EQ(wt2.size() ,n);
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(wt2[j], (typename TypeParam::value_type)text[j]);
        }
        delete [] text;
    }
}


TYPED_TEST(WtAsciiTest, DeleteTest)
{
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
