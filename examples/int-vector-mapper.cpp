#include <sdsl/int_vector_mapper.hpp>
#include <string>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 1) {
        cout << "Usage: " << argv[0] << endl;
        cout << "(1) Writes an int_vector sequentially to a file" << endl;
        cout << "(2) Streams the content from file" << endl;
        cout << "(3) Remove the file" << endl;
        return 1;
    }
    string tmp_file = "tmp_file.sdsl";
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

    // (2) memory map the content of tmp_file
    {
        int_vector_mapper<> ivm(tmp_file);
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

        if(ivm != stdv) {
        	std::cerr << "ERROR: std::vector CMP failed.";
        }
        if(ivm != iv) {
        	std::cerr << "ERROR: iv CMP failed.";
        }
        if(ivm != ivf) {
        	std::cerr << "ERROR: ivf CMP failed.";
        }
    }

    // (3) remove the file as the mapper does not do that 
    {
    	sdsl::remove(tmp_file);
    }

    {
    	auto tmp_buf = temp_file_buffer<64>::create();
    	for(const auto& val : stdv) {
    		tmp_buf.push_back(val);
    	}
        if(tmp_buf != stdv) {
        	std::cerr << "ERROR: tmp_buf CMP failed." << std::endl;
        }

        // tmp buf file is deleted automatically
    }

    return 0;
}
