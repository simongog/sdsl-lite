#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " Creates a CST and CSA for a byte file and visualizes the memory utilization during construction." << endl;
        return 1;
    }

    memory_monitor::start();

    cst_sct3<> cst;
    auto start = timer::now();
    construct(cst, argv[1], 1);
    auto stop = timer::now();
    cout << "construction cst time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;

    memory_monitor::stop();

    std::cout << "peak usage = " << memory_monitor::peak() / (1024*1024) << " MB" << std::endl;

    std::ofstream cstofs("cst-construction.html");
    cout << "writing memory usage visualization to cst-construction.html\n";
    memory_monitor::write_memory_log<HTML_FORMAT>(cstofs);
    cstofs.close();
    util::clear(cst);

    memory_monitor::start();

    csa_wt<> csa;
    start = timer::now();
    construct(csa, argv[1], 1);
    stop = timer::now();
    cout << "construction csa time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;

    memory_monitor::stop();
    std::ofstream csaofs("csa-construction.html");
    cout << "writing memory usage visualization to csa-construction.html\n";
    memory_monitor::write_memory_log<HTML_FORMAT>(csaofs);
    csaofs.close();

}
