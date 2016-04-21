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

struct bv_tag {}; // bitvector tag
struct iv_tag {}; // int_vector tag

struct csa_tag {}; // compressed suffix array (CSAs) tag
struct cst_tag {}; // compressed suffix tree (CST) tag
struct wt_tag {};  // wavelet tree tag

struct psi_tag {}; // tag for CSAs based on the psi function
struct lf_tag {}; // tag for CSAs based on the LF function

struct csa_member_tag {}; // tag for text, bwt, LF, \Psi members of CSA

struct lcp_tag {};
struct lcp_plain_tag {};
struct lcp_permuted_tag {};
struct lcp_tree_compressed_tag {};
struct lcp_tree_and_lf_compressed_tag {};

struct alphabet_tag {};
struct byte_alphabet_tag { static const uint8_t WIDTH=8; };
struct int_alphabet_tag { static const uint8_t WIDTH=0; };

struct sa_sampling_tag {};
struct isa_sampling_tag {};


template<class t_T, class t_r = void>
struct enable_if_type {
    typedef t_r type;
};

template<class t_idx, class t_enable = void>
struct index_tag {
    typedef t_enable type;
};

template<class t_idx>
struct index_tag<t_idx, typename enable_if_type<typename t_idx::index_category>::type> {
    using type = typename t_idx::index_category;
};

template<class t_sampling, class t_enable = void>
struct sampling_tag {
    typedef t_enable type;
};

template<class t_sampling>
struct sampling_tag<t_sampling, typename enable_if_type<typename t_sampling::sampling_category>::type> {
    using type = typename t_sampling::sampling_category;
};

template<class t_enc_vec, class t_enable = void>
struct is_enc_vec {
    static const bool value = false;
};

template<class t_enc_vec>
struct is_enc_vec<t_enc_vec, typename enable_if_type<typename t_enc_vec::enc_vec_type>::type> {
    static const bool value = true;
};

template<class t_alphabet, class t_enable = void>
struct is_alphabet {
    static const bool value = false;
};

template<class t_alphabet>
struct is_alphabet<t_alphabet, typename enable_if_type<typename t_alphabet::alphabet_category>::type> {
    static const bool value = true;
};

} // end namespace sdsl

#endif
