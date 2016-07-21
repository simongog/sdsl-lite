/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope thatk2_ it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file k2_treap_algorithm.hpp
    \brief k2_treap_algorithm.hpp contains k^2-treap algorithms.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_K2_TREAP_ALGORITHM
#define INCLUDED_SDSL_K2_TREAP_ALGORITHM

#include "vectors.hpp"
#include "bits.hpp"
#include "k2_tree_helper.hpp"
#include "k2_tree.hpp"
#include <tuple>
#include <algorithm>
#include <iterator>
#include <climits>
#include <vector>
#include <complex>
#include <queue>
#include <array>
#include <jmorecfg.h>

//! Namespace for the succinct data structure library.
namespace sdsl {

    namespace k2_treap_ns {

//! Check if point x is contained in the rectangle (p1,p2)
/*! \param p Point.
 *  \param Lower left corner of the rectangle.
 *  \param Upper right corner of the rectangle.
 */
        bool
        contained(const point_type p, const point_type &p1, const point_type &p2) {
            return real(p) >= real(p1) and real(p) <= real(p2) and
                   imag(p) >= imag(p1) and imag(p) <= imag(p2);
        }

//! Check if the rectangle of node v is contained in the rectangle (p1,p2)
        template<uint8_t t_k>
        bool
        contained(const point_type &p1, const point_type &p2, const node_type &v) {
//    uint64_t d = (1ULL << v.t)-1;
//    uint64_t d = (1ULL << v.t)-1;
            uint64_t d = precomp<t_k>::exp(v.t) - 1;
            return real(p1) <= real(v.p) and real(p2) >= real(v.p) + d and
                   imag(p1) <= imag(v.p) and imag(p2) >= imag(v.p) + d;
        }

//! Check if rectangle (p1,p2) and the area of node v overlap
        template<uint8_t t_k>
        bool
        overlap(const point_type &p1, const point_type &p2, const node_type &v) {
//    uint64_t d = (1ULL << v.t)-1;
            uint64_t d = precomp<t_k>::exp(v.t) - 1;
            return real(p1) <= real(v.p) + d and real(p2) >= real(v.p) and
                   imag(p1) <= imag(v.p) + d and imag(p2) >= imag(v.p);
        }
    }
// forward declaration
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    class k2_tree;


//! Specialized version of method ,,construct'' for k2_treaps.
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    void
    construct(k2_tree<t_k, t_bv, t_rank> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree<t_k, t_bv, t_rank> tmp(buf_x, buf_y, false);
        tmp.swap(idx);
    }

    //! Specialized version of method ,,construct'' for k2_treaps.
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank>
    void
    construct_bottom_up(k2_tree<t_k, t_bv, t_rank> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree<t_k, t_bv, t_rank> tmp(buf_x, buf_y, true);
        tmp.swap(idx);
    }

//! Specialized version of method ,,construct_im'' for k2_treaps.
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank,
            typename t_vector>
    void
    construct_im(k2_tree<t_k, t_bv, t_rank> &idx, t_vector data) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree<t_k, t_bv, t_rank> tmp(data, tmp_prefix,false);
        tmp.swap(idx);
    }

    //! Specialized version of method ,,construct_im'' for k2_treaps.
    template<uint8_t t_k,
            typename t_bv,
            typename t_rank,
            typename t_vector>
    void
    construct_im_bottom_up(k2_tree<t_k, t_bv, t_rank> &idx, t_vector data) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree<t_k, t_bv, t_rank> tmp(data, tmp_prefix,true);
        tmp.swap(idx);
    }
}
#endif
