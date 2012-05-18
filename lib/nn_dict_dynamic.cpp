#include "sdsl/nn_dict_dynamic.hpp"
#include "sdsl/util.hpp"

namespace sdsl
{
namespace util
{
void set_zero_bits(nn_dict_dynamic& nn)
{
    util::set_zero_bits(nn.m_tree);
}
} // end util
} // end sdsl
