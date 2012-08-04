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
 *    - rank_ith_symbol(i)
 */

#include "wt.hpp"
#include "wt_int.hpp"
#include "wt_huff.hpp"
#include "wt_rlmn.hpp"
#include "wt_rlg.hpp"
#include "wt_rlg8.hpp"

#endif
