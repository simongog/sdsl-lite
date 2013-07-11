#include <iostream>
#include <string>
#include <sdsl/util.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/louds_tree.hpp>

using namespace std;
using namespace sdsl;

typedef cst_sct3<>       cst_t;
typedef cst_t::node_type node_t;


void print_tree(const louds_tree<>& tree, const louds_tree<>::node_type& v, int depth, bit_vector& visited)
{
    typedef louds_tree<>::node_type louds_node_t;
    for (int i=0; i < depth; ++i) cout << " ";
    cout << v << "  tree.id(v) = " << tree.id(v) << endl;
    visited[tree.id(v)] = 1;
    for (uint64_t i = 1; i <= tree.degree(v); ++i) {
        louds_node_t child = tree.child(v, i);
        print_tree(tree, child, depth+1, visited);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: "<<argv[0]<< " file [max_print=80]" << std::endl;
        cout << " (1) Builds the CST of file" << std::endl;
        cout << " (2) Builds a LOUDS tree " << std::endl;
        cout << " (3) Prints information about the tree, if file size < 59." << std::endl;
        return 1;
    }
    string file = argv[1];
    cout << file << endl;
    cst_t cst;
    construct(cst, file, 1);
    uint64_t max_print = 80;
    if (argc > 2) {
        max_print = stoull(argv[2]);
    }

    typedef cst_bfs_iterator<cst_t> iterator;
    iterator begin = iterator(&cst, cst.root());
    iterator end   = iterator(&cst, cst.root(), true, true);

    louds_tree<> louds(cst, begin, end);
    if (cst.size() <= max_print+1) {
        cout << "LOUDS = " << louds.bv << endl;
        bit_vector visited(louds.nodes(), 0);
        print_tree(louds, louds.root(), 0, visited);
        cout << "visited = " << visited << endl;
    }
}
