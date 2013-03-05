#include <sdsl/suffixarrays.hpp>
#include <sdsl/lcp_construct.hpp>
#include <sdsl/bwt_construct.hpp>
#include "sdsl/config.hpp" // for CMAKE_SOURCE_DIR
#include "gtest/gtest.h"
#include <vector>
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef std::map<std::string, void (*)(sdsl::cache_config&)> tMSFP;

// The fixture for testing class int_vector.
class LcpConstructTest : public ::testing::Test
{
        protected:

                LcpConstructTest():checkprefix("CHECK_") {
                        // You can do set-up work for each test here.
						
                }

                virtual ~LcpConstructTest() {
                        // You can do clean-up work that doesn't throw exceptions here.
                }

                // If the constructor and destructor are not enough for setting up
                // and cleaning up each test, you can define the following methods:
                virtual void SetUp() {
                        // Code here will be called immediately after the constructor (right
                        // before each test).
                        std::string test_cases_dir = std::string(SDSL_XSTR(CMAKE_SOURCE_DIR)) + "/test/test_cases";
                        test_cases.push_back(sdsl::cache_config(false, test_cases_dir+"/crafted/", "example01.txt"));
                        test_cases.push_back(sdsl::cache_config(false, test_cases_dir+"/crafted/", "100a.txt"));
                        test_cases.push_back(sdsl::cache_config(false, test_cases_dir+"/crafted/", "empty.txt"));
                        // test_cases.push_back(sdsl::cache_config(false, test_cases_dir+"/crafted/", "abc_abc_abc.txt"));
                        // test_cases.push_back(sdsl::cache_config(false, test_cases_dir+"/crafted/", "abc_abc_abc2.txt")); //TODO check the problem with this one
                        test_cases.push_back(sdsl::cache_config(false, test_cases_dir+"/small/", "faust.txt"));
                        test_cases.push_back(sdsl::cache_config(false, test_cases_dir+"/small/", "zarathustra.txt"));

// 			lcp_function["construct_lcp_semi_extern_PHI"] = &sdsl::construct_lcp_semi_extern_PHI; // TODO: Handle empty test case
                        // lcp_function["construct_lcp_simple2_9n"] = &sdsl::construct_lcp_simple2_9n;
                        // lcp_function["construct_lcp_go"] = &sdsl::construct_lcp_go;
                        // lcp_function["construct_lcp_goPHI"] = &sdsl::construct_lcp_goPHI;
                        // lcp_function["construct_lcp_go2"] = &sdsl::construct_lcp_go2;
                        lcp_function["construct_lcp_bwt_based"] = &sdsl::construct_lcp_bwt_based;
                        lcp_function["construct_lcp_bwt_based2"] = &sdsl::construct_lcp_bwt_based2;
                        lcp_function["construct_lcp_PHI"] = &sdsl::construct_lcp_PHI<8>;


                        sdsl::tMSS tmp_file_map;
                        for (size_t i=0; i< this->test_cases.size(); ++i) {
							const char* file = (test_cases[i].dir+"/"+test_cases[i].id).c_str();
							uint8_t num_bytes = 1;
							{
								// Prepare Input
								sdsl::int_vector<8> text;
								ASSERT_EQ(true, sdsl::util::load_vector_from_file(text, file, num_bytes));
								ASSERT_EQ(true, sdsl::contains_no_zero_symbol(text, file) );
								sdsl::append_zero_symbol(text);
								ASSERT_EQ(true, sdsl::util::store_to_cache(text, sdsl::constants::KEY_TEXT, test_cases[i]));
								// Construct SA
								sdsl::int_vector<> sa(text.size(), 0, sdsl::bit_magic::l1BP(text.size())+1);
								sdsl::algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
								ASSERT_EQ(true, sdsl::util::store_to_cache(sa, sdsl::constants::KEY_SA, test_cases[i]));
							}
							{// Construct BWT
								sdsl::construct_bwt<8>(test_cases[i]);
							}
							{// Construct LCP
								sdsl::construct_lcp_kasai<8>(test_cases[i]);
								std::rename(sdsl::util::cache_file_name(sdsl::constants::KEY_LCP, test_cases[i]).c_str(), 
									        sdsl::util::cache_file_name((checkprefix+sdsl::constants::KEY_LCP).c_str(), test_cases[i]).c_str());
								test_cases[i].file_map.erase(sdsl::constants::KEY_LCP);
							}
                        }
                }

                virtual void TearDown() {
                        // Code here will be called immediately after each test (right
                        // before the destructor).
						for (size_t i=0; i< this->test_cases.size(); ++i) {
							std::remove(sdsl::util::cache_file_name((checkprefix+sdsl::constants::KEY_LCP).c_str(), test_cases[i]).c_str());
							sdsl::util::delete_all_files(test_cases[i].file_map);
						}
                }

                // Objects declared here can be used by all tests in the test case for Foo.
                std::vector<sdsl::cache_config> test_cases;
                tMSFP lcp_function;
				std::string checkprefix;
};

TEST_F(LcpConstructTest, construct_lcp)
{
        for (tMSFP::const_iterator it = this->lcp_function.begin(), end = this->lcp_function.end(); it != end; ++it) {
                for (size_t i=0; i< this->test_cases.size(); ++i) {
			 			std::cout << (it->first) << " on test file " << sdsl::util::cache_file_name(sdsl::constants::KEY_TEXT, this->test_cases[i]) << std::endl;

                        // Prepare LCP-Array construction
						sdsl::cache_config config_new = this->test_cases[i];

                        // Construct LCP-Array
                        (it->second)(this->test_cases[i]);

                        // Check LCP-Array
                        sdsl::int_vector<> lcp_check, lcp;
						std::string lcp_check_filename = sdsl::util::cache_file_name((checkprefix+sdsl::constants::KEY_LCP).c_str(), this->test_cases[i]);
						std::string lcp_filename = sdsl::util::cache_file_name((sdsl::constants::KEY_LCP), this->test_cases[i]);
                        ASSERT_EQ(true, sdsl::util::load_from_file(lcp_check, lcp_check_filename.c_str()))
                        << (it->first) << " on test file " << this->test_cases[i].id << " could not load reference lcp array";
                        ASSERT_EQ(true, sdsl::util::load_from_file(lcp, lcp_filename.c_str()))
                        << (it->first) << " on test file " << this->test_cases[i].id << " could not load created lcp array";
                        ASSERT_EQ(lcp_check.size(), lcp.size())
                        << (it->first) << " on test file " << this->test_cases[i].id << " lcp array size differ";
                        for (size_type j=0; j<lcp_check.size() and j<lcp.size(); ++j) {
                                ASSERT_EQ(lcp_check[j], lcp[j])
                                << (it->first) << " on test file '" << this->test_cases[i].id << "' value differ:"
                                << " lcp_check[" << j << "]=" << lcp_check[j] << "!=" << lcp[j] << "=lcp["<< j << "]";
                        }

                        // Clean up LCP-Array
                        std::remove(sdsl::util::cache_file_name(sdsl::constants::KEY_LCP, config_new).c_str());
                }
        }
}

}  // namespace

int main(int argc, char** argv)
{
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
}
