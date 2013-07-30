/* sdsl - succinct data structures library
    Copyright (C) 2010 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file algorithms_for_compressed_suffix_arrays.hpp
    \brief algorithms_for_compressed_suffix_arrays.hpp contains algorithms for compressed suffix arrays.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS
#define INCLUDED_SDSL_ALGORITHMS_FOR_COMPRESSED_SUFFIX_ARRAYS

#include "int_vector.hpp" // for bit_vector
#include "util.hpp"
#include <stack> // for calculate_supercartesian_tree_bp

namespace sdsl
{

namespace algorithm
{

template<class Csa, uint8_t int_width>
void set_isa_samples(int_vector_buffer<int_width>& sa_buf, typename Csa::isa_sample_type& isa_sample)
{
    typedef typename Csa::size_type size_type;
    size_type  n = sa_buf.size();

    isa_sample.width(bits::hi(n)+1);
    if (n >= 1) { // so n+Csa::isa_sample_dens >= 2
        isa_sample.resize((n-1+Csa::isa_sample_dens-1)/Csa::isa_sample_dens + 1);
    }
    util::set_to_value(isa_sample, 0);

    for (size_type i=0; i < n; ++i) {
        size_type sa = sa_buf[i];
        if ((sa % Csa::isa_sample_dens) == 0) {
            isa_sample[sa/Csa::isa_sample_dens] = i;
        } else if (sa+1 == n) {
            isa_sample[(sa+Csa::isa_sample_dens-1)/Csa::isa_sample_dens] = i;
        }
    }
}

}// end namespace algorithm

}// end namespace sdsl

#endif

