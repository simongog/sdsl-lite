/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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
/*! \file bp_support.hpp
    \brief bp_support.hpp contains several classed which support find_open, find_close, enclose and rr-enclose queries.
    \author Simon Gog
*/

#ifndef INCLUDED_SDSL_BP_SUPPORT
#define INCLUDED_SDSL_BP_SUPPORT

/** \defgroup bps Balanced Parentheses Supports (BPS)
 * This group contains data structures which supports a sdsl::bit_vector with the following methods:
 *   - find_open
 *   - find_close
 *   - enclose
 *   - double_enclose
 *   - rank
 *   - select
 *   - excess
 *   - rr_enclose
 */

#include "bp_support_g.hpp"
#include "bp_support_gg.hpp"
#include "bp_support_sada.hpp"

#endif
