/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file suffix_arrays.hpp
    \brief suffix_arrays.hpp contains generic classes for different suffix array classes.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SUFFIX_ARRAYS
#define INCLUDED_SDSL_SUFFIX_ARRAYS

#include "sdsl_concepts.hpp"

/** \defgroup csa Compressed Suffix Arrays (CSA) */

#include "csa_bitcompressed.hpp"
#include "csa_wt.hpp"
#include "csa_sada.hpp"
#include "wavelet_trees.hpp"
#include "construct.hpp"
#include "suffix_array_algorithm.hpp"

namespace sdsl
{

//! Typedef for convenient usage of std integer alphabet strategy
template<class t_wt              = wt_int<>,
         uint32_t t_dens         = 32,
         uint32_t t_inv_dens     = 64,
         class t_sa_sample_strat = sa_order_sa_sampling<>,
         class t_isa             = int_vector<>
         >
using csa_wt_int = csa_wt<t_wt, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, int_alphabet<>>;

template<class t_enc_vec         = enc_vector<>,          // Vector type used to store the Psi-function
         uint32_t t_dens         = 32,                    // Sample density for suffix array (SA) values
         uint32_t t_inv_dens     = 64,                    // Sample density for inverse suffix array (ISA) values
         class t_sa_sample_strat = sa_order_sa_sampling<>,// Policy class for the SA sampling. Alternative text_order_sa_sampling.
         class t_isa             = int_vector<>           // Container for the ISA samples.
         >
using csa_sada_int = csa_sada<t_enc_vec, t_dens, t_inv_dens, t_sa_sample_strat, t_isa, int_alphabet<>>;

}

#endif
