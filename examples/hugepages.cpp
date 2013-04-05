#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

template<class tCsa>
void do_something(const tCsa& csa)
{
    util::stop_watch sw;
    uint64_t sum=0;
    sw.start();
    for (size_t i=0; i<csa.size() and i<10000000; ++i) {
        sum+=csa.psi(i);
    }
    sw.stop();
    cout << sw.real_time() << endl;
    cout<<"sum="<<sum<<endl;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " (1) Creates a CST for a byte file. " << endl;
        cout << " (2) Runs a benchmark with enabled/disabled 1GB=hugepages." << endl;
        return 1;
    }
    csa_wt<> csa;
    construct(csa, argv[1], 1);
    do_something(csa); // before it is mapped
    if (mm::map_hp()) {
        cout << "Now the memory is mapped to hugepages " << endl;
        do_something(csa); // while it is mapped
        mm::unmap_hp();
    } else {
        cout << "Not able to map the memory to hugepages" << endl;
    }
    do_something(csa); // after it is unmapped
}
