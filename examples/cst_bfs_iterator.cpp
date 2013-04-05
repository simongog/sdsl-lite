#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <string>

using namespace std;
using namespace sdsl;

typedef cst_sct3<> cst_t;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "usage: "<<argv[0]<< " file" << std::endl;
        return 1;
    }

    cst_t cst;
    construct(cst, argv[1], 1);

    typedef cst_bfs_iterator<cst_t> iterator;
    iterator begin = iterator(&cst, cst.root());
    iterator end   = iterator(&cst, cst.root(), true, true);

    for (iterator it = begin; it != end; ++it) {
        std::cout << cst.depth(*it) << "-[" << cst.lb(*it) << "," << cst.rb(*it) << "]" << std::endl;
    }

}
