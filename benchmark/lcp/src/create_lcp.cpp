#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_lcp.hpp>
#include <string>
#include <chrono>
#include <iomanip>

using namespace sdsl;
using namespace std;
using namespace std::chrono;

#define S(x) #x
#define SX(x) S(x)

int main(int argc, char** argv)
{
    memory_monitor::start();
    string dir = ".";
    string id = "tmp";
    cache_config config(false, dir, id);

    register_cache_file(conf::KEY_TEXT, config);
    register_cache_file(conf::KEY_SA, config);
    register_cache_file(conf::KEY_BWT, config);

    auto start = high_resolution_clock::now();
    LCP_TYPE(config);
    auto stop = high_resolution_clock::now();
    memory_monitor::stop();
    cout << std::fixed;
    cout << "# " SX(LCPID) "_TIME = " << std::setprecision(2) << duration_cast<milliseconds>(stop-start).count()/(double)1000 << endl;
    cout << "# " SX(LCPID) "_MMPEAK = "<< memory_monitor::peak() << endl;

    return 0;
}
