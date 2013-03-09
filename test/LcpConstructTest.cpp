#include <sdsl/suffixarrays.hpp>
#include <sdsl/lcp_construct.hpp>
#include <sdsl/bwt_construct.hpp>
#include <sdsl/testutils.hpp>
#include "sdsl/config.hpp" // for CMAKE_SOURCE_DIR
#include "gtest/gtest.h"
#include <vector>
#include <string>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef std::map<std::string, void (*)(sdsl::cache_config&)> tMSFP;

// The fixture for testing class int_vector.
class LcpConstructTest : public ::testing::Test {
	protected:
	LcpConstructTest():checkprefix("CHECK_") { }

	virtual ~LcpConstructTest() { }

	// If the constructor and destructor are not enough for setting up
	// and cleaning up each test, you can define the following methods:
	virtual void SetUp() {
		sw.start();
		sw.stop();
		print_R_line(sw.get_abs_real_time(), "Kasai", "red");
		std::string prefix		= std::string(SDSL_XSTR(CMAKE_SOURCE_DIR))+"/test";
		std::string config_file = prefix + "/LcpConstructTest.config";
		std::string tc_prefix	= prefix + "/test_cases";
		std::vector<std::string> paths = sdsl::paths_from_config_file(config_file, tc_prefix.c_str());

		for ( size_t i = 0; i < paths.size(); ++i ) {
			std::string dirname = sdsl::util::dirname(paths[i]);
			std::string basename = sdsl::util::basename(paths[i]);
			test_cases.push_back( sdsl::cache_config(false, dirname, basename) );
		}
		lcp_function["construct_lcp_bwt_based"] = &sdsl::construct_lcp_bwt_based;
		lcp_function["construct_lcp_bwt_based2"] = &sdsl::construct_lcp_bwt_based2;
		lcp_function["construct_lcp_PHI"] = &sdsl::construct_lcp_PHI<8>;

		std::cout << std::fixed << std::setprecision(0);
		for (size_t i=0; i< this->test_cases.size(); ++i) {
			uint8_t num_bytes = 1;
			std::cout << "prepare " << sdsl::util::cache_file_name((checkprefix+sdsl::constants::KEY_TEXT), this->test_cases[i]) << std::endl;
			{
				// Prepare Input
				sw.start();
				std::string file = test_cases[i].dir+"/"+test_cases[i].id;
				sdsl::int_vector<8> text;
				ASSERT_EQ(true, sdsl::util::load_vector_from_file(text, file, num_bytes));
				ASSERT_EQ(true, sdsl::contains_no_zero_symbol(text, file) );
				sdsl::append_zero_symbol(text);
				ASSERT_EQ(true, sdsl::util::store_to_cache(text, sdsl::constants::KEY_TEXT, test_cases[i]));
				// Construct SA
				sdsl::int_vector<> sa(text.size(), 0, sdsl::bit_magic::l1BP(text.size())+1);
				sdsl::algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
				ASSERT_EQ(true, sdsl::util::store_to_cache(sa, sdsl::constants::KEY_SA, test_cases[i]));
				sw.stop();
				std::cout << std::setw(10) << sw.get_real_time() << "ms for sa construction" << std::endl;
			}
			{// Construct BWT
				sw.start();
				sdsl::construct_bwt<8>(test_cases[i]);
				sw.stop();
				std::cout << std::setw(10) << sw.get_real_time() << "ms for bwt construction" << std::endl;
			}
			{// Construct LCP
				sw.start();
				sw.stop();
				sw.start();
				print_R_line(sw.get_abs_real_time(), "begin"+test_cases[i].id, "gray");
				sdsl::construct_lcp_kasai<8>(test_cases[i]);
				std::rename(sdsl::util::cache_file_name(sdsl::constants::KEY_LCP, test_cases[i]).c_str(), 
							sdsl::util::cache_file_name(checkprefix+sdsl::constants::KEY_LCP, test_cases[i]).c_str());
				sw.stop();
				print_R_line(sw.get_abs_real_time(), "end  "+test_cases[i].id, "black");
				test_cases[i].file_map.erase(sdsl::constants::KEY_LCP);
				std::cout << std::setw(10) << sw.get_real_time() << "ms for lcp construction" << std::endl;
			}
		}
		print_R_line(sw.get_abs_real_time(), "Kasai", "red");
	}

	virtual void TearDown() {
		for (size_t i=0; i< this->test_cases.size(); ++i) {
			std::remove(sdsl::util::cache_file_name(checkprefix+sdsl::constants::KEY_LCP, test_cases[i]).c_str());
			sdsl::util::delete_all_files(test_cases[i].file_map);
		}
	}

	void print_R_line(uint64_t diff, std::string s, std::string color="black") {
		// std::cout << "#\t" << diff << "\t" << s << "\t" << color << std::endl;
	}
	std::vector<sdsl::cache_config> test_cases;
	tMSFP lcp_function;
	std::string checkprefix;
	sdsl::stop_watch sw;
};

TEST_F(LcpConstructTest, construct_lcp)
{
	for (tMSFP::const_iterator it = this->lcp_function.begin(), end = this->lcp_function.end(); it != end; ++it) {
			std::cout << "test" << (it->first) << std::endl;
		print_R_line(sw.get_abs_real_time(), (it->first), "red");
		for (size_t i=0; i< this->test_cases.size(); ++i) {
			// std::cout << (it->first) << " on test file " << sdsl::util::cache_file_name(sdsl::constants::KEY_TEXT, this->test_cases[i]) << std::endl;

			// Prepare LCP array construction
			sdsl::cache_config config_new = this->test_cases[i];

			// Construct LCP array
			this->sw.start();
			this->sw.stop();
			print_R_line(sw.get_abs_real_time(), "beg"+this->test_cases[i].id, "gray");
			this->sw.start();
			(it->second)(this->test_cases[i]);
			this->sw.stop();
			print_R_line(sw.get_abs_real_time(), "end"+this->test_cases[i].id, "black");
			std::cout << std::setw(10) << this->sw.get_real_time() << "ms for " << sdsl::util::cache_file_name(sdsl::constants::KEY_TEXT, this->test_cases[i]) << std::endl;

			// Check LCP array
			sdsl::int_vector<> lcp_check, lcp;
			std::string lcp_check_filename = sdsl::util::cache_file_name(checkprefix+sdsl::constants::KEY_LCP, this->test_cases[i]);
			std::string lcp_filename = sdsl::util::cache_file_name((sdsl::constants::KEY_LCP), this->test_cases[i]);
			ASSERT_EQ(true, sdsl::util::load_from_file(lcp_check, lcp_check_filename))
			<< (it->first) << " on test file " << this->test_cases[i].id << " could not load reference lcp array";
			ASSERT_EQ(true, sdsl::util::load_from_file(lcp, lcp_filename))
			<< (it->first) << " on test file " << this->test_cases[i].id << " could not load created lcp array";
			ASSERT_EQ(lcp_check.size(), lcp.size())
			<< (it->first) << " on test file " << this->test_cases[i].id << " lcp array size differ";
			for (size_type j=0; j<lcp_check.size() and j<lcp.size(); ++j) {
				ASSERT_EQ(lcp_check[j], lcp[j])
				<< (it->first) << " on test file '" << this->test_cases[i].id << "' value differ:"
				<< " lcp_check[" << j << "]=" << lcp_check[j] << "!=" << lcp[j] << "=lcp["<< j << "]";
			}
			// Clean up LCP array
			std::remove(sdsl::util::cache_file_name(sdsl::constants::KEY_LCP, config_new).c_str());
		}
		print_R_line(sw.get_abs_real_time(), (it->first), "red");
	}
}

}  // namespace

int main(int argc, char** argv)
{
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
}

