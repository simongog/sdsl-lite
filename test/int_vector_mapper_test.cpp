#include "sdsl/int_vector_mapper.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <random>
#include <algorithm>

namespace
{

typedef sdsl::int_vector<>::size_type size_type;
typedef sdsl::int_vector<>::value_type value_type;

std::string temp_dir;

// The fixture for testing class int_vector.
class int_vector_mapper_test : public ::testing::Test
{
    protected:

        int_vector_mapper_test() {}

        virtual ~int_vector_mapper_test() {}

        virtual void SetUp()
        {
            std::mt19937_64 rng;
            {
                std::uniform_int_distribution<uint64_t> distribution(1, 100000);
                auto dice = bind(distribution, rng);
                for (size_type i=0; i<10; ++i) {
                    vec_sizes.push_back(dice());
                }
            }
        }

        virtual void TearDown() {}

        std::vector<size_type> vec_sizes = {1,64,65,127,128}; // different sizes for the vectors
};

TEST_F(int_vector_mapper_test, iterator)
{
    static_assert(std::is_move_constructible<sdsl::int_vector_mapper<>>::value, "Type is not move constructible");
    static_assert(std::is_move_assignable<sdsl::int_vector_mapper<>>::value, "Type is not move assignable");

    // test plain
    for (const auto& size : vec_sizes) {
        std::vector<uint64_t> vec(size);
        sdsl::util::set_to_id(vec);
        {
            std::ofstream ofs(temp_dir+"/int_vector_mapper_itrtest",std::ios::binary | std::ios::trunc | std::ios::out);
            sdsl::serialize_vector(vec,ofs);
        }
        {
            const sdsl::int_vector_mapper<64,std::ios_base::in> ivm(temp_dir+"/int_vector_mapper_itrtest",true);
            ASSERT_EQ(size,ivm.size());
            ASSERT_TRUE(std::equal(ivm.begin(),ivm.end(),vec.begin()));
            ASSERT_EQ(size,(size_t)std::distance(ivm.begin(),ivm.end()));
        }
        {
            const sdsl::int_vector_mapper<64,std::ios_base::in> ivm(temp_dir+"/int_vector_mapper_itrtest",true);
            auto itr = ivm.end()-1;
            for (size_t i=0; i<size; i++) {
                ASSERT_EQ(*itr,vec[size-i-1]);
                --itr;
            }
        }
        sdsl::remove(temp_dir+"/int_vector_mapper_itrtest");
    }

    // test fixed width
    for (const auto& size : vec_sizes) {
        sdsl::int_vector<25> vec(size);
        sdsl::util::set_to_id(vec);
        store_to_file(vec,temp_dir+"/int_vector_mapper_itrtest");
        {
            const sdsl::int_vector_mapper<25,std::ios_base::in> ivm(temp_dir+"/int_vector_mapper_itrtest");
            ASSERT_EQ(size,ivm.size());
            ASSERT_TRUE(std::equal(ivm.begin(),ivm.end(),vec.begin()));
            ASSERT_EQ(size,(size_t)std::distance(ivm.begin(),ivm.end()));
        }
        {
            const sdsl::int_vector_mapper<25,std::ios_base::in> ivm(temp_dir+"/int_vector_mapper_itrtest");
            auto itr = ivm.end()-1;
            for (size_t i=0; i<size; i++) {
                ASSERT_EQ(*itr,vec[size-i-1]);
                --itr;
            }
        }
        sdsl::remove(temp_dir+"/int_vector_mapper_itrtest");
    }

    // test variable width
    for (const auto& size : vec_sizes) {
        sdsl::int_vector<> vec(size);
        sdsl::util::set_to_id(vec);
        sdsl::util::bit_compress(vec);
        store_to_file(vec,temp_dir+"/int_vector_mapper_itrtest");
        {
            const sdsl::int_vector_mapper<0,std::ios_base::in> ivm(temp_dir+"/int_vector_mapper_itrtest");
            ASSERT_EQ(size,ivm.size());
            ASSERT_EQ(vec.width(),ivm.width());
            ASSERT_TRUE(std::equal(ivm.begin(),ivm.end(),vec.begin()));
            ASSERT_EQ(size,(size_t)std::distance(ivm.begin(),ivm.end()));
        }
        {
            sdsl::int_vector_mapper<> ivm(temp_dir+"/int_vector_mapper_itrtest");
            auto itr = ivm.end()-1;
            for (size_t i=0; i<size; i++) {
                ASSERT_EQ(*itr,vec[size-i-1]);
                --itr;
            }
        }
        sdsl::remove(temp_dir+"/int_vector_mapper_itrtest");
    }
}

TEST_F(int_vector_mapper_test, push_back)
{
    // test plain
    for (const auto& size : vec_sizes) {
        std::vector<uint64_t> vec(size);
        sdsl::util::set_to_id(vec);
        {
            std::ofstream ofs(temp_dir+"/int_vector_mapper_push_backtest",std::ios::binary | std::ios::trunc | std::ios::out);
            sdsl::serialize_vector(vec,ofs);
        }
        {
            sdsl::int_vector_mapper<64> ivm(temp_dir+"/int_vector_mapper_push_backtest",true);
            ASSERT_EQ(size,ivm.size());
            for (size_t i=0; i<size; i++) {
                vec.push_back(vec.size());
                ivm.push_back(ivm.size());
            }
        }
        {
            sdsl::int_vector_mapper<64> ivm(temp_dir+"/int_vector_mapper_push_backtest",true);
            ASSERT_EQ(vec.size(),ivm.size());
            ASSERT_TRUE(std::equal(ivm.begin(),ivm.end(),vec.begin()));
            ASSERT_EQ(vec.size(),(size_t)std::distance(ivm.begin(),ivm.end()));
        }
        sdsl::remove(temp_dir+"/int_vector_mapper_itrtest");
    }

    // test fixed width
    for (const auto& size : vec_sizes) {
        sdsl::int_vector<31> vec(size);
        std::vector<uint64_t> stdvec(size);
        sdsl::util::set_to_id(vec);
        sdsl::util::set_to_id(stdvec);
        store_to_file(vec,temp_dir+"/int_vector_mapper_push_backtest");
        {
            sdsl::int_vector_mapper<31> ivm(temp_dir+"/int_vector_mapper_push_backtest");
            ASSERT_EQ(size,ivm.size());
            for (size_t i=0; i<size; i++) {
                stdvec.push_back(stdvec.size());
                ivm.push_back(ivm.size());
            }
        }
        {
            sdsl::int_vector_mapper<31> ivm(temp_dir+"/int_vector_mapper_push_backtest");
            ASSERT_EQ(stdvec.size(),ivm.size());
            ASSERT_TRUE(std::equal(ivm.begin(),ivm.end(),stdvec.begin()));
            ASSERT_EQ(stdvec.size(),(size_t)std::distance(ivm.begin(),ivm.end()));
        }
        sdsl::remove(temp_dir+"/int_vector_mapper_push_backtest");
    }

    // test variable width
    for (const auto& size : vec_sizes) {
        sdsl::int_vector<> vec(size);
        std::vector<uint64_t> stdvec(size);
        sdsl::util::set_to_id(vec);
        sdsl::util::set_to_id(stdvec);
        sdsl::util::bit_compress(vec);
        store_to_file(vec,temp_dir+"/int_vector_mapper_push_backtest");
        {
            sdsl::int_vector_mapper<> ivm(temp_dir+"/int_vector_mapper_push_backtest");
            ASSERT_EQ(size,ivm.size());
            for (size_t i=0; i<size; i++) {
                stdvec.push_back(i);
                ivm.push_back(i);
            }
        }
        {
            sdsl::int_vector_mapper<> ivm(temp_dir+"/int_vector_mapper_push_backtest");
            ASSERT_EQ(stdvec.size(),ivm.size());
            ASSERT_TRUE(std::equal(ivm.begin(),ivm.end(),stdvec.begin()));
            ASSERT_EQ(stdvec.size(),(size_t)std::distance(ivm.begin(),ivm.end()));
        }
        sdsl::remove(temp_dir+"/int_vector_mapper_push_backtest");
    }
}

TEST_F(int_vector_mapper_test, bit_compress)
{
    for (const auto& size : vec_sizes) {
        sdsl::int_vector<> vec(size);
        sdsl::util::set_to_id(vec);
        store_to_file(vec,temp_dir+"/int_vector_mapper_bit_compress_test_uncompressed");
        sdsl::util::bit_compress(vec);
        store_to_file(vec,temp_dir+"/int_vector_mapper_bit_compress_test");
        {
            sdsl::int_vector_mapper<> ivmc(temp_dir+"/int_vector_mapper_bit_compress_test");
            sdsl::int_vector_mapper<> ivmu(temp_dir+"/int_vector_mapper_bit_compress_test_uncompressed");
            ASSERT_TRUE(std::equal(ivmc.begin(),ivmc.end(),ivmu.begin()));
            ASSERT_TRUE(std::equal(ivmc.begin(),ivmc.end(),vec.begin()));
        }
        {
            sdsl::int_vector_mapper<> ivmu(temp_dir+"/int_vector_mapper_bit_compress_test_uncompressed");
            ASSERT_TRUE(std::equal(ivmu.begin(),ivmu.end(),vec.begin()));
            sdsl::util::bit_compress(ivmu);
            ASSERT_TRUE(std::equal(ivmu.begin(),ivmu.end(),vec.begin()));
        }
        {
            sdsl::int_vector_mapper<> ivmc(temp_dir+"/int_vector_mapper_bit_compress_test");
            sdsl::int_vector_mapper<> ivmu(temp_dir+"/int_vector_mapper_bit_compress_test_uncompressed");
            ASSERT_EQ(ivmc.size(),ivmu.size());
            ASSERT_EQ(ivmc.width(),ivmu.width());
            ASSERT_TRUE(std::equal(ivmc.begin(),ivmc.end(),ivmu.begin()));
        }
        sdsl::remove(temp_dir+"/int_vector_mapper_bit_compress_test_uncompressed");
        sdsl::remove(temp_dir+"/int_vector_mapper_bit_compress_test");
    }
}

TEST_F(int_vector_mapper_test, bitvector_mapping)
{
    for (const auto& size : vec_sizes) {
        sdsl::bit_vector bv(size);
        sdsl::util::set_random_bits(bv,4711);
        store_to_file(bv,temp_dir+"/bit_vector_mapper_test");
        {
            // load/store test
            const sdsl::bit_vector_mapper<std::ios_base::in> bvm(temp_dir+"/bit_vector_mapper_test");
            ASSERT_EQ(bvm.size(),bv.size());
            ASSERT_EQ(bvm.width(),bv.width());
            ASSERT_TRUE(std::equal(bvm.begin(),bvm.end(),bv.begin()));
            ASSERT_EQ(sdsl::util::cnt_one_bits(bv),sdsl::util::cnt_one_bits(bvm));
        }
        {
            // flip test
            sdsl::bit_vector_mapper<> bvm(temp_dir+"/bit_vector_mapper_test");
            bvm.flip();
            bv.flip();
            ASSERT_TRUE(std::equal(bvm.begin(),bvm.end(),bv.begin()));
            ASSERT_EQ(sdsl::util::cnt_one_bits(bv),sdsl::util::cnt_one_bits(bvm));
        }
        {
            // load/store after flip
            const sdsl::bit_vector_mapper<std::ios_base::in> bvm(temp_dir+"/bit_vector_mapper_test");
            ASSERT_TRUE(std::equal(bvm.begin(),bvm.end(),bv.begin()));
            ASSERT_EQ(sdsl::util::cnt_one_bits(bv),sdsl::util::cnt_one_bits(bvm));
        }
        sdsl::remove(temp_dir+"/bit_vector_mapper_test");
    }
}

TEST_F(int_vector_mapper_test, temp_buffer_test)
{
    for (const auto& size : vec_sizes) {
        sdsl::int_vector<> vec(size);
        sdsl::util::set_to_id(vec);
        std::string tmp_file_name;
        {
            auto tmp_buf = sdsl::temp_file_buffer<31>::create();
            tmp_file_name = tmp_buf.file_name();
            ASSERT_EQ(tmp_buf.width(),(uint8_t)31);
            ASSERT_EQ(tmp_buf.size(),(size_t)0);
            ASSERT_TRUE(tmp_buf.empty());
            for (const auto& val : vec) {
                tmp_buf.push_back(val);
            }
            ASSERT_EQ(tmp_buf.size(),vec.size());
            ASSERT_TRUE(std::equal(tmp_buf.begin(),tmp_buf.end(),vec.begin()));
        }
        // check that the file is gone
        std::ifstream cfs(tmp_file_name);
        ASSERT_FALSE(cfs.is_open());
    }
}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 2) {
        // LCOV_EXCL_START
        std::cout << "Usage: " << argv[0] << " tmp_dir" << std::endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    temp_dir = argv[1];
    return RUN_ALL_TESTS();
}
