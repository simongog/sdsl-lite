#include <iostream>
#include <fstream>
#include <sdsl/construct.hpp>

using namespace sdsl;
using namespace std;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char** argv)
{
    // Check Parameters
    if (argc!=2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " Computes the SA from file space efficiently." << endl;
        cout << " Result is stored in `sa_<file>.sdsl`." << endl;
        return 1;
    }

    cache_config config(false, ".", util::basename(argv[1]));

    int_vector<8> text;
    load_vector_from_file(text, argv[1], 1);
    append_zero_symbol(text);
    auto n = text.size();
    store_to_cache(text, conf::KEY_TEXT, config);
    util::clear(text);

    cout << "Calculate suffix array ... " << flush;
    construct_config::byte_algo_sa = SE_SAIS; // or LIBDIVSUFSORT for less space-efficient but faster construction
    memory_monitor::start();
    auto start = timer::now();
    construct_sa<8>(config);
    auto stop = timer::now();
    memory_monitor::stop();
    cout << "done." << endl;
    cout << "Construction needed:" << endl;
    cout << duration_cast<seconds>(stop-start).count() << " seconds." << endl;
    cout << memory_monitor::peak() << " bytes." << endl;
    cout << (1.0*memory_monitor::peak())/n << " bytes per input byte." << endl;
    sdsl::remove(cache_file_name(conf::KEY_TEXT, config));
    return 0;
}
