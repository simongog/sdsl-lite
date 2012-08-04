#include <iostream>
#include <string>
#include <sdsl/util.hpp>
#include <sdsl/suffixtrees.hpp>
#include <sdsl/louds_tree.hpp>

using namespace std;
using namespace sdsl;

template<class tTree, class tNode>
void print_tree(const tTree& tree, const tNode& v, int depth, bit_vector& visited)
{
    typedef typename tTree::size_type size_type;
    if (tree.nodes() < 60) {
        for (int i=0; i<depth; ++i) cout << " ";
        cout << v << "  tree.id(v) = "<<tree.id(v)<<endl;
        visited[tree.id(v)] = 1;
    }
    for (size_type i = 1; i <= tree.degree(v); ++i) {
        tNode child = tree.child(v, i);
        print_tree(tree, child, depth+1, visited);
        tNode parent = tree.parent(child);
        if (parent != v) {
            cout << "ERROR: tree.parent("<<child<<")="<< parent << "!=" << v <<endl;
        }
    }
}

template<class tCst>
void test(string file)
{
    std::cout << file << std::endl;
    tCst cst;
//	util::verbose = true;
    construct_cst(file, cst);

    typedef cst_bfs_iterator<tCst> iterator;
    iterator begin = iterator(&cst, cst.root());
    iterator end   = iterator(&cst, cst.root(), true, true);

    louds_tree<> louds(cst, begin, end);
    if (louds.nodes() < 60) {
        cout << "LOUDS = " << louds.bv << endl;
        bit_vector visited(louds.nodes(), 0);
        print_tree(louds, louds.root(), 0, visited);
        cout << "visited = " << visited << endl;
    }
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "usage: "<<argv[0]<< " file_name" << std::endl;
    } else {
        test<cst_sct3<> >(argv[1]);
    }
}
