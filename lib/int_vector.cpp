#include "sdsl/int_vector.hpp"

namespace sdsl
{


template<>
void int_vector<0>::width(const uint8_t new_int_width)
{
    if (0 < new_int_width and new_int_width <= 64)
        m_width = new_int_width;
    else
        m_width = 64;
}

template<>
void bit_vector::flip()
{
    if (!empty()) {
        for (uint64_t i=0; i<(capacity()>>6); ++i) {
            m_data[i] = ~m_data[i];
        }
    }
}

}
