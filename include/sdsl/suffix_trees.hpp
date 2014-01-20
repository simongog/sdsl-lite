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
/*! \file suffix_trees.hpp
    \brief suffix_trees.hpp contains generic classes for different suffix tree classes.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SUFFIX_TREES
#define INCLUDED_SDSL_SUFFIX_TREES

#include "sdsl_concepts.hpp"
#include "suffix_arrays.hpp"
#include "suffix_tree_algorithm.hpp"
#include "construct.hpp"
#include "util.hpp"
#include <iostream>
#include <cmath>
#include <set>

using std::cout;
using std::endl;

namespace sdsl
{

// Gets ISA[SA[idx]+d]
// d = depth of the character 0 = first position
template<class t_csa>
typename t_csa::size_type get_char_pos(typename t_csa::size_type idx, typename t_csa::size_type d, const t_csa& csa)
{
    if (d == 0)
        return idx;
    // if we have to apply \f$\LF\f$ or \f$\Phi\f$ more
    // than 2*d times to calc csa(csa[idx]+d), we opt to
    // apply \f$ \Phi \f$ d times
    if (csa.sa_sample_dens + csa.isa_sample_dens > 2*d+2) {
        for (typename t_csa::size_type i=0; i < d; ++i)
            idx = csa.psi[idx];
        return idx;
    }
    return csa.isa[csa[idx] + d];
}

}

/** \defgroup cst Compressed Suffix Trees (CST)
 *   This group contains data structures for compressed suffix trees. The following methods are supported:
 *    - root()
 *    - child(v,c)
 *    - select_child(v)
 *    - select_leaf(i)
 *    - parent(v)
 *    - sl(v)
 *    - lca(v,w)
 *    - ..
 */

#include "suffix_tree_helper.hpp"
#include "cst_sct3.hpp"
#include "cst_sada.hpp"

#include "csa_bitcompressed.hpp"
#include "int_vector.hpp"

#include <iostream>
#include <string>

#endif
