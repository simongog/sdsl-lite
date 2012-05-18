#ifndef INCLUDED_SDSL_WT_HELPER
#define INCLUDED_SDSL_WT_HELPER

#include "int_vector.hpp"

namespace sdsl
{

const uint16_t _undef_node = 65535;

//! Count for each character in [0..255] the number of occurences in rac[0..size-1]
/*!
 * \return C An array of size 256, which contains for each character the number of occurences in rac[0..size-1]
 */
template<class size_type_class, class size_type>
void calculate_character_occurences(int_vector_file_buffer<8, size_type_class>& text, const size_type size, size_type* C)
{
    text.reset();
    for (size_type i=0, r_sum=0, r = text.load_next_block(); r_sum < size;) {
        for (; i < r_sum+r; ++i) {
            ++C[text[i-r_sum]];
        }
        r_sum += r; r = text.load_next_block();
    }
}

template<class size_type>
void calculate_effective_alphabet_size(const size_type* C, size_type& sigma)
{
    sigma = 0; // initialize with 0
    for (size_type i=0; i < 256; ++i) // for all possible symbols
        if (C[i] > 0)				  // increase sigma, if it
            ++sigma;				 // exists in the text
}

} // end namespace sdsl
#endif
