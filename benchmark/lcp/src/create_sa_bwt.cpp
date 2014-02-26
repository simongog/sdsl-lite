#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_sa.hpp>
#include <sdsl/construct_bwt.hpp>
#include <string>
#include <chrono>
#include <iostream>

using namespace sdsl;
using namespace std;
using namespace std::chrono;

typedef bit_vector::size_type size_type;

//argv[1] = test file
int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout<<"Usage: input_file" << std::endl;
    }
    memory_monitor::start();
    string file = argv[1];
    uint8_t num_bytes = 1; // Byte Alphabet
    string dir = ".";
    string id = "tmp";
    cache_config config(false, dir, id);

    //load text
    auto start = high_resolution_clock::now();
    {
        int_vector<8> text;
        load_vector_from_file(text, file, num_bytes);
        if (contains_no_zero_symbol(text, file)) {
            append_zero_symbol(text);
            store_to_cache(text, conf::KEY_TEXT, config);
        }
        register_cache_file(conf::KEY_TEXT, config);
    }
    auto stop = high_resolution_clock::now();
    memory_monitor::stop();
    cout << "# TXT_TIME = " << duration_cast<milliseconds>(stop-start).count()/(double)1000 << endl;
    cout << "# TXT_MMPEAK = " << memory_monitor::peak() << endl;

    //construct sa
    memory_monitor::start();
    start = high_resolution_clock::now();
    {
        construct_sa<8>(config);
        register_cache_file(conf::KEY_SA, config);
    }
    stop = high_resolution_clock::now();
    memory_monitor::stop();
    cout << "# SA_TIME = " << duration_cast<milliseconds>(stop-start).count()/(double)1000 << endl;
    cout << "# SA_MMPEAK = " << memory_monitor::peak() << endl;

    //construct bwt
    memory_monitor::start();
    start = high_resolution_clock::now();
    {
        construct_bwt<8>(config);
        register_cache_file(conf::KEY_BWT, config);
    }
    stop = high_resolution_clock::now();
    memory_monitor::stop();
    cout << "# BWT_TIME = " << duration_cast<milliseconds>(stop-start).count()/(double)1000 <<endl;
    cout << "# BWT_MMPEAK = "<< memory_monitor::peak() << endl;

    return 0;
}

