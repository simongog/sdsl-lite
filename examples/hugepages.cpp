#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

template<class t_cst>
void do_something(const t_cst& cst)
{
    uint64_t sum=0;
    auto start = timer::now();
    for (size_t i=0; i<cst.csa.size() and i<50000000; ++i) {
        sum+=cst.csa.lf[i];
    }
    auto stop = timer::now();
    cout << "runtime in seconds: " << duration_cast<seconds>(stop-start).count() << endl;
    cout << "sum="<<sum<<endl;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file <use_hugepages>" << endl;
        cout << " (1) Creates a CST for a byte file. " << endl;
        cout << " (2) Runs a benchmark with enabled/disabled hugepages." << endl;
        return 1;
    }

    if (argc==3) {
        // memory_manager::use_hugepages(500ULL*1024ULL*1024ULL);
        // use all available hugepages if nothing is specified
        memory_manager::use_hugepages();
    }

    cst_sct3<> cst;
    auto start = timer::now();
    construct(cst, argv[1], 1);
    auto stop = timer::now();
    cout << "construction time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;
    do_something(cst); // before it is mapped
}
