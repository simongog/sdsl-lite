#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <string>

using namespace std;
using namespace sdsl;

typedef cst_sct3<> cst_t;
typedef cst_sada<> csts_t;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "usage: "<<argv[0]<< " file" << std::endl;
        return 1;
    }

    cst_t cst;
    construct(cst, argv[1], 1);

    auto root = cst.root();

    for (auto& child: cst.children(root)) {
        std::cout << "sct3 id = " << cst.id(child) << std::endl;
    }

    csts_t csts;
    construct(csts, argv[1], 1);
    auto roots = csts.root();
    for (auto child: csts.children(roots)) {
        std::cout << "sada id = " << csts.id(child) << std::endl;
    }
}
