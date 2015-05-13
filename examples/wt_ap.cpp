#include <sdsl/wavelet_trees.hpp>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        return 1;
    }
    wt_ap<> wt;
    memory_monitor::start();
    construct(wt, argv[1], 1);
    memory_monitor::stop();
    {
        std::ofstream csaofs("wt-construction.html");
        cout << "writing memory usage visualization to csa-construction.html\n";
        memory_monitor::write_memory_log<HTML_FORMAT>(csaofs);
    }
    {
        std::ofstream out("wt-space.html");
        write_structure<HTML_FORMAT>(wt,out);
    }
}
