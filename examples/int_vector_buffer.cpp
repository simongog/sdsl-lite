#include <sdsl/int_vector_buffer.hpp>
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

    // (1) Writes an int_vector sequentially to the file tmp_file
    {
        // create an int_vector_buffer
        int_vector_buffer<> ivb(tmp_file,      // filename
                                std::ios::out, // we do not want to open an existing file, but create a new one
                                1024*1024,     // use a buffer of about 1MB
                                64,            // use 64bit for each integer
                                false);        // use int_vector format

        // write sequentially random values to disk
        for (uint64_t i=0; i<size; ++i) {
            ivb[i] = rng();
        }
    }

    // (2) Streams the content of tmp_file
    {
        // open file with an int_vector_buffer.
        // Default is: open existing file, use 1MB Buffer, 64bit for each integer and int_vector format
        int_vector_buffer<> ivb(tmp_file);
        if (!ivb.is_open()) {
            std::cerr << "ERROR: ivb could not be opend with file " << tmp_file << size << std::endl;
            return 1;
        }
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
        ivb.close(true); // close buffer and (3) remove the file
    }

    return 0;
}
