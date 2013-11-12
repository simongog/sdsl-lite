#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_sa.hpp>
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
    string file = argv[1];
    uint8_t num_bytes = 1; // Byte Alphabet
    string dir = ".";
    string id = "tmp";
    cache_config config(false, dir, id);

    //load text
    auto start = high_resolution_clock::now();
    int_vector<8> text;
    load_vector_from_file(text, file, num_bytes);
    if (contains_no_zero_symbol(text, file)) {
        append_zero_symbol(text);
        store_to_cache(text, conf::KEY_TEXT, config);
    }
    register_cache_file(conf::KEY_TEXT, config);
    auto stop = high_resolution_clock::now();
    cout <<"# TXT_TIME = " << duration_cast<milliseconds>(stop-start).count() <<endl;

    //construct sa

    start = high_resolution_clock::now();
    construct_sa<8>(config);
    register_cache_file(conf::KEY_SA, config);
    stop = high_resolution_clock::now();
    cout <<"# SA_TIME = "<<duration_cast<seconds>(stop-start).count() << endl;
    int i = system("echo -n \"# SA_VMPEAK = \";grep VmPeak /proc/self/status | grep -o '[0-9]*'");
    i = system("echo -n \"# SA_VMHWM = \";grep VmPeak /proc/self/status | grep -o '[0-9]*'");
    return i;
}
