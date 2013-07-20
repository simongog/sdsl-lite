#include "sdsl/suffix_arrays.hpp"
#include "gtest/gtest.h"
#include <cstdlib>
#include <vector>
#include <string>

namespace
{

using namespace sdsl;
using namespace std;

typedef int_vector<>::size_type size_type;
tMSS    test_case_file_map;
string  test_file;
uint8_t num_bytes;
string  temp_file;
string  temp_dir;
bool in_memory;


template<class T>
class CsaIntTest : public ::testing::Test { };


using testing::Types;

typedef Types<  csa_wt<wt_int<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, int_alphabet<> >,
        csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, int_alphabet<> >,
        csa_bitcompressed<int_alphabet<> >,
        csa_wt<wt_int<rrr_vector<63> >, 8, 8, sa_order_sa_sampling<>, int_vector<>, int_alphabet<> >,
        csa_wt<wt_int<>, 32, 32, text_order_sa_sampling<>, int_vector<>, int_alphabet<> >,
        csa_sada<enc_vector<>, 32, 32, text_order_sa_sampling<>, int_vector<>, int_alphabet<> >
        > Implementations;

TYPED_TEST_CASE(CsaIntTest, Implementations);

TYPED_TEST(CsaIntTest, CreateAndStoreTest)
{
    TypeParam csa;
    cache_config config(false, temp_dir, util::basename(test_file));
    construct(csa, test_file, config, num_bytes);
    test_case_file_map = config.file_map;
    bool success = store_to_file(csa, temp_file);
    ASSERT_EQ(true, success);
}

//! Test access methods
TYPED_TEST(CsaIntTest, Sigma)
{
    TypeParam csa;
    ASSERT_EQ(true, load_from_file(csa, temp_file));
    int_vector<> text;
    load_vector_from_file(text, test_file, num_bytes);
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

//! Test suffix array access methods
TYPED_TEST(CsaIntTest, SaAccess)
{
    TypeParam csa;
    ASSERT_EQ(true, load_from_file(csa, temp_file));
    int_vector<> sa;
    load_from_file(sa, test_case_file_map[constants::KEY_SA]);
    size_type n = sa.size();
    ASSERT_EQ(n, csa.size());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ(sa[j], csa[j])<<" j="<<j;
    }
}


//! Test inverse suffix access methods
TYPED_TEST(CsaIntTest, IsaAccess)
{
    TypeParam csa;
    ASSERT_EQ(true, load_from_file(csa, temp_file));
    int_vector<> isa;
    size_type n = 0;
    {
        int_vector<> sa;
        load_from_file(sa, test_case_file_map[constants::KEY_SA]);
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

//! Test Burrows-Wheeler access methods
TYPED_TEST(CsaIntTest, BwtAccess)
{
    if (test_case_file_map.end() != test_case_file_map.find(constants::KEY_BWT_INT)) {
        TypeParam csa;
        ASSERT_EQ(true, load_from_file(csa, temp_file));
        int_vector<> bwt;
        load_from_file(bwt, test_case_file_map[constants::KEY_BWT_INT]);
        size_type n = bwt.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(bwt[j], csa.bwt[j])<<" j="<<j;
        }
    }
}

//! Test text access methods
TYPED_TEST(CsaIntTest, TextAccess)
{
    if (test_case_file_map.end() != test_case_file_map.find(constants::KEY_TEXT_INT)) {
        TypeParam csa;
        ASSERT_EQ(true, load_from_file(csa, temp_file));
        int_vector<> text;
        load_from_file(text, test_case_file_map[constants::KEY_TEXT_INT]);
        size_type n = text.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(text[j], csa.text[j])<<" j="<<j;
        }
    }
}

//! Test Psi access methods
TYPED_TEST(CsaIntTest, PsiAccess)
{
    if (test_case_file_map.end() != test_case_file_map.find(constants::KEY_PSI)) {
        TypeParam csa;
        ASSERT_EQ(true, load_from_file(csa, temp_file));
        int_vector<> psi;
        load_from_file(psi, test_case_file_map[constants::KEY_PSI]);
        size_type n = psi.size();
        ASSERT_EQ(n, csa.size());
        for (size_type j=0; j<n; ++j) {
            ASSERT_EQ(psi[j], csa.psi[j])<<" j="<<j;
        }
    }
}

//! Test if Psi[LF[i]]=i
TYPED_TEST(CsaIntTest, PsiLFAccess)
{
    TypeParam csa;
    ASSERT_EQ(true, load_from_file(csa, temp_file));
    for (size_type j=0; j<csa.size(); ++j) {
        size_type lf = csa.psi(j);
        ASSERT_TRUE(lf >= 0);
        ASSERT_TRUE(lf < csa.size());
        ASSERT_EQ(j, csa.psi[lf])<<" j="<<j;
    }
}


//! Test access after swap
TYPED_TEST(CsaIntTest, SwapTest)
{
    TypeParam csa1;
    ASSERT_EQ(true, load_from_file(csa1, temp_file));
    TypeParam csa2;
    csa1.swap(csa2);
    int_vector<> sa;
    load_from_file(sa, test_case_file_map[constants::KEY_SA]);
    size_type n = sa.size();
    ASSERT_EQ(n, csa2.size());
    for (size_type j=0; j<n; ++j) {
        ASSERT_EQ((typename TypeParam::value_type)sa[j], csa2[j]);
    }
}

TYPED_TEST(CsaIntTest, DeleteTest)
{
    std::remove(temp_file.c_str());
    util::delete_all_files(test_case_file_map);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " test_file num_bytes temp_file tmp_dir" << endl;
        cout << " (1) Generates a CSA out of test_file; stores it in temp_file." << endl;
        cout << "     Temporary files (SA/BWT/TEXT) are stored in tmp_dir." << endl;
        cout << "     num_bytes specifies who many bytes make a symbol in the"<< endl;
        cout << "     input sequence" << endl;
        cout << "     If `in-memory` is specified, the in-memory construction is tested." << endl;
        cout << " (2) Performs tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
    }
    test_file = argv[1];
    num_bytes = atoi(argv[2]);
    temp_file = argv[3];
    temp_dir  = argv[4];
    in_memory    = argc > 5;
    if (in_memory) {
        temp_dir = "@";
        int_vector<> data;
        load_vector_from_file(data, test_file, num_bytes);
        test_file = ram_file_name(test_file);
        switch (num_bytes) {
            case 0: store_to_file(data, test_file); break;
            case 1: store_to_plain_array<uint8_t>(data, test_file); break;
            case 2: store_to_plain_array<uint16_t>(data, test_file); break;
            case 3: store_to_plain_array<uint32_t>(data, test_file); break;
            case 4: store_to_plain_array<uint64_t>(data, test_file); break;
        }
        temp_file = ram_file_name(temp_file);
    }
    return RUN_ALL_TESTS();
}
