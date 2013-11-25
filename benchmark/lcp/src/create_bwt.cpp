#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_bwt.hpp>
#include <string>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace sdsl;
using namespace std;
using namespace std::chrono;

int main(int argc, char** argv)
{
    memory_monitor::start();
    string dir = ".";
    string id = "tmp";
    cache_config config(false, dir, id);

    register_cache_file(conf::KEY_TEXT, config);
    register_cache_file(conf::KEY_SA, config);

    auto start = high_resolution_clock::now();
    construct_bwt<8>(config);
    register_cache_file(conf::KEY_BWT, config);
    auto stop = high_resolution_clock::now();
    memory_monitor::stop();
    cout << std::fixed;
    cout << "# BWT_TIME = " << std::setprecision(2) << duration_cast<milliseconds>(stop-start).count()/(double)1000 <<endl;
    cout << "# BWT_MMPEAK = "<< memory_monitor::peak() << endl;

    return 0;
}
