#include <iostream>
#include <string>
#include <sdsl/util.hpp>
#include <sdsl/suffix_trees.hpp>

using namespace std;
using namespace sdsl;


template<class tCst>
void test(const char* file) {
    std::cout << file << std::endl;
    tCst cst;
//	util::verbose = true;
    construct(cst, file, 1);

    typedef cst_bfs_iterator<tCst> iterator;
    iterator begin = iterator(&cst, cst.root());
    iterator end   = iterator(&cst, cst.root(), true, true);

    for (iterator it = begin; it != end; ++it) {
        std::cout << cst.depth(*it) << "-[" << cst.lb(*it) << "," << cst.rb(*it) << "]" << std::endl;
    }
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "usage: "<<argv[0]<< " file" << std::endl;
    } else {
        test<cst_sct3<> >(argv[1]);
    }
}
