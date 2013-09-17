#include <sdsl/suffix_arrays.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/construct_bwt.hpp>
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <map>

using namespace sdsl;
using namespace std;

namespace
{

string test_file, temp_dir, test_id;
typedef map<string, void (*)(cache_config&)> tMSFP;// map <name, lcp method>

// The fixture for testing class int_vector.
class LcpConstructTest : public ::testing::Test
{
    protected:
        LcpConstructTest():CHECK_KEY("CHECK_"+string(conf::KEY_LCP)) { }

        virtual ~LcpConstructTest() { }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {
            test_config = cache_config(false, temp_dir, test_id);
            lcp_function["bwt_based"] = &construct_lcp_bwt_based;
            lcp_function["bwt_based2"] = &construct_lcp_bwt_based2;
            lcp_function["PHI"] = &construct_lcp_PHI<8>;
            lcp_function["semi_extern_PHI"] = &construct_lcp_semi_extern_PHI;
            lcp_function["go"] = &construct_lcp_go;
            lcp_function["goPHI"] = &construct_lcp_goPHI;

            uint8_t num_bytes = 1;
            {
                // Prepare Input
                int_vector<8> text;
                ASSERT_TRUE(load_vector_from_file(text, test_file, num_bytes));
                ASSERT_TRUE(contains_no_zero_symbol(text, test_file));
                append_zero_symbol(text);
                ASSERT_TRUE(store_to_cache(text, conf::KEY_TEXT, test_config));
                // Construct SA
                int_vector<> sa(text.size(), 0, bits::hi(text.size())+1);
                algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
                ASSERT_TRUE(store_to_cache(sa, conf::KEY_SA, test_config));
            }
            {
                // Construct BWT
                construct_bwt<8>(test_config);
            }
            {
                // Construct LCP
                construct_lcp_kasai<8>(test_config);
                std::rename(cache_file_name(conf::KEY_LCP, test_config).c_str(),
                            cache_file_name(CHECK_KEY, test_config).c_str());
                test_config.file_map.erase(conf::KEY_LCP);
            }
        }

        virtual void TearDown() {
            sdsl::remove(cache_file_name(CHECK_KEY, test_config));
        }

        cache_config test_config;
        tMSFP lcp_function;
        string CHECK_KEY;
};

TEST_F(LcpConstructTest, construct_lcp)
{
    for (tMSFP::const_iterator it = this->lcp_function.begin(), end = this->lcp_function.end(); it != end; ++it) {
        string info = "construct_lcp_" + (it->first) + " on test file " + test_file;
        // Construct LCP array
        (it->second)(this->test_config);
        // Check LCP array
        int_vector<> lcp_check, lcp;
        string lcp_check_file = cache_file_name(CHECK_KEY, this->test_config);
        string lcp_file = cache_file_name(conf::KEY_LCP, this->test_config);
        ASSERT_TRUE(load_from_file(lcp_check, lcp_check_file))
                << info << " could not load reference lcp array";
        ASSERT_TRUE(load_from_file(lcp, lcp_file))
                << info << " could not load created lcp array";
        ASSERT_EQ(lcp_check.size(), lcp.size())
                << info << " lcp array size differ";
        for (uint64_t j=0; j<lcp.size(); ++j) {
            ASSERT_EQ(lcp_check[j], lcp[j])
                    << info << " value differ:" << " lcp_check[" << j << "]="
                    << lcp_check[j] << "!=" << lcp[j] << "=lcp["<< j << "]";
        }
        // Clean up LCP array
        sdsl::remove(cache_file_name(conf::KEY_LCP, this->test_config));
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 4) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " test_file tmp_dir ID" << endl;
        cout << " (1) Generates the SA, BWT and LCP; arrays are stored in tmp_dir." << endl;
        cout << "     File contain ID as substring." << endl;
        cout << " (2) Generates LCP with other algorithm and checks the result." << endl;
        cout << " (3) Deletes all generated files." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    test_file = argv[1];
    temp_dir  = argv[2];
    test_id   = argv[3];
    return RUN_ALL_TESTS();
}
