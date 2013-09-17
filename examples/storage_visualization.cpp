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
        sum+=csa.lf[i];
    }
    auto stop = timer::now();
    cout << "runtime in s: " << duration_cast<seconds>(stop-start).count() << endl;
    cout <<"sum="<<sum<<endl;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " Creates a CST for a byte file and visualizes the space used by the data structure." << endl;
        return 1;
    }

    cst_sct3<> cst;
    auto start = timer::now();
    construct(cst, argv[1], 1);
    auto stop = timer::now();
    cout << "construction cst time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;

    std::ofstream cstofs("cst-space-usage.html");
    cout << "writing storage visualization to cst-space-usage.html\n";
    write_structure<HTML_FORMAT>(cst,cstofs);
    cstofs.close();
    util::clear(cst);

}
