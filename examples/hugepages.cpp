#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

template<class tCsa>
void do_something(const tCsa& csa)
{
    uint64_t sum=0;
    auto start = timer::now();
    for (size_t i=0; i<csa.size() and i<10000000; ++i) {
        sum+=csa.psi(i);
    }
    auto stop = timer::now();
    cout << "runtime in ms: " << duration_cast<microseconds>(stop-start).count() << endl;
    cout <<"sum="<<sum<<endl;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " (1) Creates a CST for a byte file. " << endl;
        cout << " (2) Runs a benchmark with enabled/disabled 1GB=hugepages." << endl;
        return 1;
    }

    if (argc==3) {
        memory_manager::use_hugepages(500*1024*1024);
    }

    csa_wt<> csa;
    auto start = timer::now();
    construct(csa, argv[1], 1);
    auto stop = timer::now();
    cout << "construction in ms: " << duration_cast<microseconds>(stop-start).count() << endl;
    do_something(csa); // before it is mapped
}
