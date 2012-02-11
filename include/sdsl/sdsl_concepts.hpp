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
/*! \file suffixarray_helper.hpp
    \brief suffixarray_helper.hpp contains some helper classes for (compressed suffix arrays)
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CONCEPTS
#define INCLUDED_SDSL_CONCEPTS

namespace sdsl{

struct csa_tag{};

struct cst_tag{};

struct lcp_plain_tag{};
struct lcp_permuted_tag{};
struct lcp_tree_compressed_tag{};
struct lcp_tree_and_lf_compressed_tag{};

} // end namespace sdsl

#endif
