#include <iostream>
#include <sdsl/wt_topk.hpp>
#include <sdsl/wt_topk_algorithm.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    wt_topk<> wtk;
    construct_im(wtk, {{0,0,2},{1,2,3},{2,1,2},{3,0,2},{4,0,1},{5,1,4},{6,0,1},{7,1,1},{8,0,8},{9,2,5}});
    wtk.print_info();
}
