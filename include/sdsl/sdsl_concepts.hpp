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
/*! \file sdsl_concepts.hpp
    \brief Contains declarations and definitions of data structure concepts.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_CONCEPTS
#define INCLUDED_SDSL_CONCEPTS

#include "uintx_t.hpp" // for uint8_t

namespace sdsl
{

struct csa_tag {}; // compressed suffix array (CSAs) tag
struct cst_tag {}; // compressed suffix tree (CST) tag
struct wt_tag {};  // wavelet tree tag

struct psi_tag {}; // tag for CSAs based on the psi function
struct lf_tag {}; // tag for CSAs based on the LF function

struct lcp_plain_tag {};
struct lcp_permuted_tag {};
struct lcp_tree_compressed_tag {};
struct lcp_tree_and_lf_compressed_tag {};


struct byte_alphabet_tag {
    static const uint8_t WIDTH=8;
};
struct int_alphabet_tag {
    static const uint8_t WIDTH=0;
};

} // end namespace sdsl


/*
  Define enable_if (which is now in C++11),
  and is_same.
  Reference:
    * http://en.cppreference.com/w/cpp/types/enable_if
    * http://stackoverflow.com/questions/4354665/function-that-takes-an-stl-iterator-over-any-container-of-a-elements-of-a-specif
*/

template <bool, typename T>
struct enable_if;

template <typename T>
struct enable_if<true, T> {
    typedef T type;
};

template <typename T, typename U>
struct is_same {
    enum {value = false};
};

template <typename T>
struct is_same<T, T> {
    enum {value = true};
};



#endif
