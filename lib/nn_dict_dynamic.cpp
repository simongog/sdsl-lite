#include "sdsl/nn_dict_dynamic.hpp"
#include "sdsl/util.hpp"

namespace sdsl
{
namespace util
{
void set_zero_bits(nn_dict_dynamic& nn)
{
    util::set_to_value(nn.m_tree, 0);
}
} // end util
} // end sdsl
