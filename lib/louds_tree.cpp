#include "sdsl/louds_tree.hpp"

namespace sdsl
{
std::ostream& operator<<(std::ostream& os, const louds_node& v)
{
    os<<"("<<v.nr<<","<<v.pos<<")";
    return os;
}
}
