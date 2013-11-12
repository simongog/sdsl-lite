#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_lcp.hpp>
#include <string>
#include <chrono>

using namespace sdsl;
using namespace std;
using namespace std::chrono;

#define S(x) #x
#define SX(x) S(x)

int main(int argc, char** argv)
{
    string dir = ".";
    string id = "tmp";
    cache_config config(false, dir, id);

    register_cache_file(conf::KEY_TEXT, config);
    register_cache_file(conf::KEY_SA, config);
    register_cache_file(conf::KEY_BWT, config);

    auto start = high_resolution_clock::now();
    LCP_TYPE(config);
    auto stop = high_resolution_clock::now();
    cout << "# " SX(LCPID) "_TIME = " << duration_cast<seconds>(stop-start).count() << endl;
    int i = system("echo -n \"# " SX(LCPID) "_VMPEAK = \";grep VmPeak /proc/self/status | grep -o '[0-9]*'");
    i = system("echo -n \"# " SX(LCPID) "_VMHWM = \";grep VmPeak /proc/self/status | grep -o '[0-9]*'");
    return i;
}
