#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

int main(int, char* argv[])
{
    cst_sct3<> cst;
    construct(cst, argv[1], 1);
    uint64_t max_depth = stoull(argv[2]);

    // use the DFS iterator to traverse `cst`
    for (auto it=cst.begin(); it!=cst.end(); ++it) {
        if (it.visit() == 1) {  // node visited the first time
            auto v = *it;       // get the node by dereferencing the iterator
            if (cst.depth(v) <= max_depth) {   // if depth node is <= max_depth
                // process node, e.g. output it in format d-[lb, rb]
                cout<<cst.depth(v)<<"-["<<cst.lb(v)<< ","<<cst.rb(v)<<"]"<<endl;
            } else { // skip the subtree otherwise
                it.skip_subtree();
            }
        }
    }
}
