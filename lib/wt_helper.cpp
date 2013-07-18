#include "sdsl/wt_helper.hpp"

namespace sdsl
{
pc_node::pc_node(uint64_t freq, uint64_t sym, uint64_t parent,
                 uint64_t child_left, uint64_t child_right):
    freq(freq), sym(sym), parent(parent)
{
    child[0] = child_left;
    child[1] = child_right;
}

pc_node& pc_node::operator=(const pc_node& v)
{
    freq     = v.freq;
    sym      = v.sym;
    parent   = v.parent;
    child[0] = v.child[0];
    child[1] = v.child[1];
    return *this;
}

}
