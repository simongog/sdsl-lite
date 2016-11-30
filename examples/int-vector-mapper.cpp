#include <sdsl/int_vector_mapper.hpp>
#include <string>
#include <iostream>

using namespace sdsl;
using namespace std;
                            

int main(int argc, char* argv[])
{
    std::cout << "main() " << std::endl;
    if (argc < 1) {
        cout << "Usage: " << argv[0] << endl;
        cout << "(1) Writes an int_vector sequentially to a file" << endl;
        cout << "(2) Streams the content from file" << endl;
        cout << "(3) Remove the file" << endl;
        return 1;
    }
    std::cout << "A " << std::endl;
    string tmp_file = "@tmp_file.sdsl";
    size_t size = 10000000;
    std::mt19937_64 rng(13);
    uint8_t width = 0;
    int_vector<> iv(size,0,64);
    int_vector<64> ivf(size,0);
    std::vector<uint64_t> stdv(size,0);

    // (1) write an int vector to disk
    {
        // write sequentially random values to disk
        for (uint64_t i=0; i<size; ++i) {
            iv[i] = rng();
            stdv[i] = iv[i];
            ivf[i] = iv[i];
        }

        util::bit_compress(iv);
        width = iv.width();
        store_to_file(iv,tmp_file);
    }
    std::cout << "store to file " << tmp_file << std::endl;

    // (2) open readonly! memory map the content of tmp_file
    {
        std::cout << "start mapper " << tmp_file << std::endl;
        const int_vector_mapper<0,std::ios_base::in> ivm(tmp_file);
        std::cout << "stop mapper " << tmp_file << std::endl;
        if (ivm.size() != size) {
            std::cerr << "ERROR: ivm.size()="<< ivm.size() << " != " << size << std::endl;
            return 1;
        }
        if (ivm.width() != width) {
            std::cerr << "ERROR: ivm.width()="<< ivm.width() << " != " << width << std::endl;
            return 1;
        }
        rng.seed(13); // To get the same values than before use the same seed
        for (uint64_t i=0; i<ivm.size(); ++i) {
            uint64_t expected_value = rng();
            if (ivm[i] != expected_value) {
                std::cerr << "ERROR: ivm["<< i << "]=" << ivm[i] << " != " << expected_value << "= expected_value" << std::endl;
                return 1;
            }
        }

        if (ivm != stdv) {
            std::cerr << "ERROR: std::vector CMP failed.";
        }
        if (ivm != iv) {
            std::cerr << "ERROR: iv CMP failed.";
        }
        if (ivm != ivf) {
            std::cerr << "ERROR: ivf CMP failed.";
        }
    }

    // (2) open read+write! memory map the content of tmp_file
    {
        int_vector_mapper<0> ivm(tmp_file);
        if (ivm.size() != size) {
            std::cerr << "ERROR: ivm.size()="<< ivm.size() << " != " << size << std::endl;
            return 1;
        }
        if (ivm.width() != width) {
            std::cerr << "ERROR: ivm.width()="<< ivm.width() << " != " << width << std::endl;
            return 1;
        }
        rng.seed(13); // To get the same values than before use the same seed
        for (uint64_t i=0; i<ivm.size(); ++i) {
            uint64_t expected_value = rng();
            if (ivm[i] != expected_value) {
                std::cerr << "ERROR: ivm["<< i << "]=" << ivm[i] << " != " << expected_value << "= expected_value" << std::endl;
                return 1;
            }
        }

        if (ivm != stdv) {
            std::cerr << "ERROR: std::vector CMP failed.";
        }
        if (ivm != iv) {
            std::cerr << "ERROR: iv CMP failed.";
        }
        if (ivm != ivf) {
            std::cerr << "ERROR: ivf CMP failed.";
        }
    }

    // (3) remove the file as the mapper does not do that if we don't specify it
    {
        sdsl::remove(tmp_file);
    }

    int_vector<> v(1000000,0,16);
    for(size_t i=0; i<v.size(); ++i){
        v[i] = i;
    }

    {
        auto tmp_buf = write_out_mapper<0>::create("@test",v.size(),v.width());
        std::cout<<"tmp_buf.size()="<<tmp_buf.size()<<" v.size()="<<v.size()<<std::endl;
        std::cout<<"tmp_buf.width()="<<(size_t)tmp_buf.width()<<" v.width()="<<(size_t)v.width()<<std::endl;
        for (size_t i=0; i<v.size(); ++i) {
            tmp_buf[i] = v[i];
            if ( tmp_buf[i] != v[i] ) {
                std::cout<<"i="<<i<<" tmp_buf[i]="<<tmp_buf[i]<<" != "<<v[i]<<std::endl;
                break;
            }
        }
        if (tmp_buf != v) {
            std::cerr << "ERROR: tmp_buf CMP failed." << std::endl;
        }
    }

    {
        std::cout<<"file_size="<<ram_fs::file_size("@test")<<std::endl;
        int_vector<> vv;
        load_from_file(vv, "@test");
        std::cout<<"v.size()="<<v.size()<<" ? "<<vv.size()<<std::endl;
        for(size_t i=0; i<v.size(); ++i){
            if (v[i] !=  vv[i]) {
                std::cout<<"i="<<i<<"v[i]="<<v[i]<<" != "<<vv[i]<<std::endl;
                break;
            }
        }
    }

    return 0;
}
