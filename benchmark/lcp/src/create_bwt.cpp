#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_bwt.hpp>
#include <string>
#include <chrono>

using namespace sdsl;
using namespace std;
using namespace std::chrono;

int main(int argc, char** argv)
{
    string dir = ".";
    string id = "tmp";
    cache_config config(false, dir, id);

    register_cache_file(conf::KEY_TEXT, config);
    register_cache_file(conf::KEY_SA, config);

    auto start = high_resolution_clock::now();
    construct_bwt<8>(config);
    register_cache_file(conf::KEY_BWT, config);
    auto stop = high_resolution_clock::now();
    cout << "# BWT_TIME = " << duration_cast<seconds>(stop-start).count() <<endl;
    int i = system("echo -n \"# BWT_VMPEAK = \";grep VmPeak /proc/self/status | grep -o '[0-9]*'");
    i = system("echo -n \"# BWT_VMHWM = \";grep VmPeak /proc/self/status | grep -o '[0-9]*'");
    return i;
}
