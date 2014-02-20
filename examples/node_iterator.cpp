#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

template<class t_cst>
void output_node(const typename t_cst::node_type& v, const t_cst& cst)
{
    cout << cst.depth(v) << "-[" << cst.lb(v) << ","
         << cst.rb(v) << "]" << endl;
}

template<class t_cst>
void run()
{
    t_cst cst;
    construct_im(cst, "ananas", 1);
    for (auto v : cst) {
        output_node(v, cst);
    }
    cout<<"--"<<endl;
    auto v = cst.select_leaf(2);
    for (auto it = cst.begin(v); it != cst.end(v); ++it) {
        output_node(*it, cst);
    }
    cout<<"--"<<endl;
    v = cst.parent(cst.select_leaf(4));
    for (auto it = cst.begin(v); it != cst.end(v); ++it) {
        output_node(*it, cst);
    }
    cout<<"---"<<endl;
}

int main()
{
    run<cst_sct3<>>();
    run<cst_sada<>>();
}
