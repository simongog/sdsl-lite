/* sdsl - succinct data structures library
    Copyright (C) 2011-2014 Simon Gog

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
/*! \file lcp_dac.hpp
    \brief lcp_dac.hpp contains an implementation of a (compressed) LCP array.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_LCP_DAC
#define INCLUDED_SDSL_LCP_DAC

#include "lcp.hpp"
#include "vectors.hpp"
#include "rank_support_v5.hpp"

namespace sdsl
{

//! A class for the compressed version of LCP information of an suffix array
/*! A dac_vector is used to compress represent the values compressed.
 *  The template parameter are forwarded to the dac_vector.
 *  \tparam t_b    Split block size.
 *  \tparam t_rank Rank structure to navigate between the different levels.
 */
template<uint8_t  t_b    = 4,
         typename t_rank = rank_support_v5<>>
using lcp_dac = lcp_vlc<dac_vector<t_b, t_rank>>;

template<typename t_bv = bit_vector, int t_default_max_levels = 64>
using lcp_dac_dp = lcp_vlc<dac_vector_dp<t_bv, t_default_max_levels>>;

} // end namespace sdsl
#endif
