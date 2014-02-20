/* sdsl - succinct data structures library
    Copyright (C) 2011 Simon Gog

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
/*! \file wavelet_trees.hpp
    \brief wavelet_trees.hpp contains wavelet tree implementations.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_WAVELET_TREES
#define INCLUDED_SDSL_WAVELET_TREES

/** \defgroup wt Wavelet Trees (WT)
 *   This group contains data structures for wavelet trees. The following methods are supported:
 *    - []-operator
 *    - rank(i, c)
 *    - select(i, c)
 *    - inverse_select(i)
 */

#include "wt_pc.hpp"
#include "wt_blcd.hpp"
#include "wt_gmr.hpp"
#include "wt_huff.hpp"
#include "wt_hutu.hpp"
#include "wt_int.hpp"
#include "wm_int.hpp"
#include "wt_rlmn.hpp"
#include "construct.hpp"
#include "wt_algorithm.hpp"

namespace sdsl
{

template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type
         >
using wt_hutu_int = wt_pc<hutu_shape,
      t_bitvector,
      t_rank,
      t_select,
      t_select_zero,
      int_tree<>>;

template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type>
using wt_huff_int = wt_pc<huff_shape,
      t_bitvector,
      t_rank,
      t_select,
      t_select_zero,
      int_tree<>>;

template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select_one  = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type>
using wt_blcd_int = wt_pc<balanced_shape,
      t_bitvector,
      t_rank,
      t_select_one,
      t_select_zero,
      int_tree<>>;
}

#endif
