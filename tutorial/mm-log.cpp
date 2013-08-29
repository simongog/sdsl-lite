#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <chrono>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    // set granularity of logging to 20 milliseconds
    memory_monitor::granularity(std::chrono::milliseconds(20));

    // generate CST
    memory_monitor::start();
    {
        cst_sct3<> cst;
        construct(cst, argv[1], 1);
        cerr<<cst.size()<<endl;
    }
    memory_monitor::stop();

    std::ofstream of("test.csv");
    memory_monitor::write_memory_log<CSV>(of);
    of.close();
}
