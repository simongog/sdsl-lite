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
#include "k2_tree.hpp"
#include "k2_tree_partitioned.hpp"
#include "k2_tree_hybrid.hpp"
#include <tuple>
#include <algorithm>
#include <iterator>
#include <climits>
#include <vector>
#include <complex>
#include <queue>
#include <array>

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
            typename t_lev,
            typename t_leaf,
            bool t_comp,
            uint8_t t_access_shortcut_size,
            typename t_rank>
    class k2_tree;

    // forward declaration
    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev,
            typename t_leaf,
            bool t_comp,
            uint8_t t_access_shortcut_size,
            typename t_rank>
    class k2_tree_hybrid;

    // forward declaration
    template<uint8_t t_k0,
            typename subk2_tree,
            bool t_comp>
    class k2_tree_partitioned;


    template<uint8_t t_k0,
            typename subk2_tree,
            bool t_comp>
    void
    construct(k2_tree_partitioned<t_k0, subk2_tree, t_comp> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree_partitioned<t_k0, subk2_tree, t_comp> tmp(buf_x, buf_y, false);
        tmp.swap(idx);
    }

    template<uint8_t t_k,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            bool t_comp=false,
            uint8_t t_access_shortcut_size=0,
            typename t_rank=typename t_lev::rank_1_type>
    void
    construct(k2_tree <t_k, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree<t_k, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> tmp(buf_x, buf_y, false);
        tmp.swap(idx);
    }

    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            bool t_comp=false,
            uint8_t t_access_shortcut_size=0,
            typename t_rank=typename t_lev::rank_1_type>
    void
    construct(k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp,t_access_shortcut_size,t_rank> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp,t_access_shortcut_size,t_rank> tmp(buf_x, buf_y, false);
        tmp.swap(idx);
    }

    template<uint8_t t_k,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type,
            bool t_comp=false,
            uint8_t t_access_shortcut_size=0>
    void
    construct_bottom_up(k2_tree<t_k, t_lev, t_leaf, t_comp,t_access_shortcut_size,t_rank> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree<t_k, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> tmp(buf_x, buf_y, true);
        tmp.swap(idx);
    }

    template<uint8_t t_k0,
            typename subk2_tree,
            bool t_comp>
    void
    construct_bottom_up(k2_tree_partitioned<t_k0, subk2_tree, t_comp> &idx, std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree_partitioned<t_k0, subk2_tree, t_comp> tmp(buf_x, buf_y, true);
        tmp.swap(idx);
    }

    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            bool t_comp=false,
            uint8_t t_access_shortcut_size=0,
            typename t_rank=typename t_lev::rank_1_type>
    void
    construct_bottom_up(k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp,t_access_shortcut_size,t_rank> &idx,
                        std::string file) {
        int_vector_buffer<> buf_x(file + ".x", std::ios::in);
        int_vector_buffer<> buf_y(file + ".y", std::ios::in);
        k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp,t_access_shortcut_size,t_rank> tmp(buf_x, buf_y, true);
        tmp.swap(idx);
    }

    template<uint8_t t_k,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type,
            typename t_vector,
            bool t_comp=false,
            uint8_t t_access_shortcut_size=0>
    void
    construct_im(k2_tree<t_k, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> &idx, t_vector data, uint64_t max_hint = 0) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree<t_k, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> tmp(tmp_prefix,false,data, max_hint);
        tmp.swap(idx);
    }

    template<uint8_t t_k0,
            typename subk2_tree,
            bool t_comp,
            typename t_vector>
    void
    construct_im(k2_tree_partitioned<t_k0, subk2_tree, t_comp> &idx, t_vector data, uint64_t max_hint = 0) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree_partitioned<t_k0, subk2_tree, t_comp> tmp(tmp_prefix,false,data, max_hint);
        tmp.swap(idx);
    }

    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type,
            typename t_vector,
            bool t_comp=false,
            uint8_t t_access_shortcut_size=0>
    void
    construct_im(k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> &idx,
                 t_vector data, uint64_t max_hint = 0) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>
                tmp(tmp_prefix,false,data, max_hint);
        tmp.swap(idx);
    }

    template<uint8_t t_k,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type,
            typename t_vector,
            bool t_comp=false,
            uint8_t t_access_shortcut_size=0>
    void
    construct_im_bottom_up(k2_tree<t_k, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> &idx, t_vector data, uint64_t max_hint = 0) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree<t_k, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> tmp(tmp_prefix,true,data, max_hint);
        tmp.swap(idx);
    }

    template<uint8_t t_k0,
            typename subk2_tree,
            bool t_comp,
            typename t_vector>
    void
    construct_im_bottom_up(k2_tree_partitioned<t_k0, subk2_tree, t_comp> &idx, t_vector data, uint64_t max_hint = 0) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree_partitioned<t_k0, subk2_tree, t_comp> tmp(tmp_prefix,true,data, max_hint);
        tmp.swap(idx);
    }

    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type,
            typename t_vector,
            uint8_t t_access_shortcut_size=0,
            bool t_comp=false>
    void
    construct_im_bottom_up(k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank> &idx,
                           t_vector data, uint64_t max_hint = 0) {
        std::string tmp_prefix = ram_file_name("k2_treap_");
        k2_tree_hybrid<t_k_l_1, t_k_l_1_size, t_k_l_2, t_k_leaves, t_lev, t_leaf, t_comp, t_access_shortcut_size, t_rank>
                tmp(tmp_prefix,true, data, max_hint);
        tmp.swap(idx);
    }
}
#endif