#ifndef INCLUDED_SDSL_WT_HELPER
#define INCLUDED_SDSL_WT_HELPER

#include "int_vector.hpp"

namespace sdsl
{

const uint16_t _undef_node = 65535;

//! Count for each character in [0..255] the number of occurrences in rac[0..size-1]
/*!
 * \return C An array of size 256, which contains for each character the number of occurrences in rac[0..size-1]
 */
void calculate_character_occurences(int_vector_file_buffer<8>& text, const int_vector_size_type size, int_vector_size_type* C);


template<class size_type, class sigma_type>
void calculate_effective_alphabet_size(const size_type* C, sigma_type& sigma)
{
    sigma = 0; // initialize with 0
    for (size_type i=0; i < 256; ++i) // for all possible symbols
        if (C[i] > 0)				  // increase sigma, if it
            ++sigma;				 // exists in the text
}

} // end namespace sdsl
#endif
