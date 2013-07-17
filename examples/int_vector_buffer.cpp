#include <sdsl/int_vector_buffer.hpp>
#include <string>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " temp_file" << endl;
        cout << "(1) Writes an int_vector sequentially to disk" << endl;
        cout << "(2) Streams the content of temp_file" << endl;
        cout << "(3) Remove temp_file" << endl;
        return 1;
    }
    string tmp_file = argv[1];
    size_t size = 10000000;
    std::mt19937_64 rng(13);

    // (1) Writes an int_vector sequentially to disk
    {
        // create an int_vector_buffer
        int_vector_buffer<> ivb(tmp_file, // filename
                                false,     // we do not want to open an existing file, but create a new one
                                1024*1024, // use a buffer of about 1MB
                                64,        // use 64bit for each integer
                                true,      // keep file when int_vector_buffer is destroyed
                                false);    // use int_vector format

        // write sequentially random values to disk
        for (uint64_t i=0; i<size; ++i) {
            ivb[i] = rng();
        }
    }

    // (2) Streams the content of temp_file
    {
        // open file with an int_vector_buffer
        int_vector_buffer<> ivb(tmp_file, // filename
                                true,           // we want to open an existing file
                                1024,           // use a buffer of about 1KB
                                64,             // use 64bit for each integer
                                false,          // do not keep file when int_vector_buffer is destroyed
                                false);         // use int_vector format
        if (ivb.size() != size) {
            std::cerr << "ERROR: ivb.size()="<< ivb.size() << " != " << size << std::endl;
            return 1;
        }
        rng.seed(13); // To get the same values than before use the same seed
        for (uint64_t i=0; i<ivb.size(); ++i) {
            uint64_t expected_value = rng();
            if (ivb[i] != expected_value) {
                std::cerr << "ERROR: ivb["<< i << "]=" << ivb[i] << " != " << expected_value << "= expected_value" << std::endl;
                return 1;
            }
        }
    } // temporary file will be removed here

    return 0;
}
