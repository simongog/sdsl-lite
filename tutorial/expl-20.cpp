#include <sdsl/suffix_arrays.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    memory_monitor::start();
    csa_wt<> csa;
    construct(csa, "english.200MB", 1);
    memory_monitor::stop();
    memory_monitor::write_memory_log<HTML_FORMAT>(cout);
}
