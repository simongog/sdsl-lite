#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    // set granularity of logging to 20 milliseconds
    mm::log_granularity(20);
    // connect cout to the logging stream
    mm::log_stream(&cout);
    // generate CST
    mm::log("begin");
    {
        cst_sct3<> cst;
        construct(cst, argv[1], 1);
        cerr<<cst.size()<<endl;
    }
    mm::log("end");
}
