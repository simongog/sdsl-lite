#include <sdsl/suffix_arrays.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/construct_sa_se.hpp>
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <map>

using namespace sdsl;
using namespace std;

namespace
{
cache_config config;
uint64_t n;

class SaConstructTest : public ::testing::Test
{
    protected:
        SaConstructTest() { }

        virtual ~SaConstructTest() { }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:
        virtual void SetUp() {}

        virtual void TearDown() {}
};

TEST_F(SaConstructTest, divsufsort)
{
    // Construct SA with divsufsort
    memory_monitor::start();
    construct_config::byte_algo_sa = LIBDIVSUFSORT;
    construct_sa<8>(config);
    memory_monitor::stop();
    // Rename sa-file
    sdsl::rename(cache_file_name(conf::KEY_SA, config), cache_file_name("check_sa", config));
    config.file_map.erase(conf::KEY_SA);
    register_cache_file("check_sa", config);
    cout << "# constructs_space = " << (1.0*memory_monitor::peak())/n << " byte per byte, =>" << memory_monitor::peak() << " bytes in total" << endl;
}

TEST_F(SaConstructTest, sesais)
{
    // Construct SA with seSAIS
    memory_monitor::start();
    construct_config::byte_algo_sa = SE_SAIS;
    construct_sa<8>(config);
    memory_monitor::stop();
    cout << "# constructs_space = " << (1.0*memory_monitor::peak())/n << " byte per byte, =>" << memory_monitor::peak() << " bytes in total" << endl;
}

TEST_F(SaConstructTest, compare)
{
    // Load both SAs
    int_vector_buffer<> sa_check(cache_file_name("check_sa", config));
    int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, config));

    // Verify
    ASSERT_EQ(sa_check.size(), sa.size()) << " suffix array size differ";
    for (uint64_t i=0; i<sa_check.size(); ++i) {
        ASSERT_EQ(sa_check[i], sa[i]) << " sa differs at position " << i;
    }

    // Remove all files
    util::delete_all_files(config.file_map);
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 4) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " test_file tmp_dir ID" << endl;
        cout << " (1) Generates the SA and stores it in tmp_dir." << endl;
        cout << "     All generated files contains ID as substring." << endl;
        cout << " (2) Checks the suffix array." << endl;
        cout << " (3) Deletes all generated files." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    string test_file = argv[1];
    string temp_dir  = argv[2];
    string test_id   = argv[3];
    config = cache_config(false, temp_dir, test_id);
    if (!cache_file_exists(conf::KEY_TEXT, config)) {
        int_vector<8> text;
        load_vector_from_file(text, test_file, 1);
        if (contains_no_zero_symbol(text, test_file)) {
            append_zero_symbol(text);
            store_to_cache(text, conf::KEY_TEXT, config);
        }
        n = text.size();
    }
    register_cache_file(conf::KEY_TEXT, config);

    return RUN_ALL_TESTS();
}
