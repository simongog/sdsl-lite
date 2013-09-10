/* sdsl - succinct data structures library
    Copyright (C) 2010-2013 Simon Gog

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
/*! \file wt_huff.hpp
    \brief wt_huff.hpp contains a class for a Huffman shaped wavelet tree
                       over byte sequences.
    \author Simon Gog and Timo Beller
*/
#ifndef INCLUDED_SDSL_WT_HUFF
#define INCLUDED_SDSL_WT_HUFF

#include "wt_pc.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

// forward declaration
struct huff_shape;

//! A Huffman-shaped wavelet tree.
/*!
 * A wavelet tree is build for a vector of characters over the byte alphabet
 * \f$\Sigma\f$. If you need a wavelet tree for a integer alphabet you should
 * use `wt_int`.
 * The wavelet tree \f$wt\f$ consists of a tree of bitvectors and provides
 * three efficient methods:
 *   - The "[]"-operator: \f$wt[i]\f$ returns the i-th symbol of vector for
 *     which the wavelet tree was build for.
 *   - The rank method: \f$wt.rank(i,c)\f$ returns the number of occurrences
 *     of symbol \f$c\f$ in the prefix [0..i-1] in the vector for which the
 *     wavelet tree was build for.
 *   - The select method: \f$wt.select(j,c)\f$ returns the index
 *     \f$i\in [0..size()-1]\f$ of the j-th occurrence of symbol \f$c\f$.
 *
 * The idea of using a Huffman shaped wavelet was first mentioned on page 17
 * of the following technical report:
 * Veli Mäkinen and Gonzalo Navarro:
 * ,,Succinct Suffix Arrays based on Run-Length Encoding.''
 * Available under: http://swp.dcc.uchile.cl/TR/2005/TR_DCC-2005-004.pdf
 *
 * \tparam t_bitvector   Underlying bitvector structure.
 * \tparam t_rank        Rank support for pattern `1` on the bitvector.
 * \tparam t_select      Select support for pattern `1` on the bitvector.
 * \tparam t_select_zero Select support for pattern `0` on the bitvector.
 * \tparam t_dfs_shape   Layout of the tree structure in memory. Set 0
 *                       for BFS layout and 1 fro DFS layout.
 *
 * \par Space complexity
 *      \f$n H_0 + 2|\Sigma|\log n\f$ bits, where \f$n\f$ is the size
 *       of the vector the wavelet tree was build for.
 *
 *  @ingroup wt
 */
template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_tree_strat  = byte_tree<> >
using wt_huff = wt_pc<huff_shape,
      t_bitvector,
      t_rank,
      t_select,
      t_select_zero,
      t_tree_strat>;

// Huffman shape for wt_pc
template<class t_wt>
struct _huff_shape {
    typedef typename t_wt::size_type         size_type;
    typedef std::pair<size_type, size_type>  tPII;    // (freq, nodenr)-pair
    typedef std::priority_queue
    <tPII, std::vector<tPII>,
    std::greater<tPII>>                      tMPQPII; // min priority queue
    enum { lex_ordered = 0 };

    template<class t_rac>
    static void
    construct_tree(t_rac& C, std::vector<pc_node>& temp_nodes) {
        tMPQPII pq;
        size_type i = 0;
        // add leaves of Huffman tree
        std::for_each(std::begin(C), std::end(C), [&](decltype(*std::begin(C)) &freq) {
            if (freq > 0) {
                pq.push(tPII(freq, temp_nodes.size()));// push (frequency, node pointer)
                // initial bv_pos with number of occurrences and bv_pos_rank
                // value with the code of the corresponding char, parent,
                // child[0], and child[1] are set to undef
                temp_nodes.emplace_back(pc_node(freq, i));
            }
            ++i;
        });
        while (pq.size() > 1) {
            tPII v1, v2;
            v1 = pq.top(); pq.pop();
            v2 = pq.top(); pq.pop();
            temp_nodes[v1.second].parent = temp_nodes.size(); // parent is new node
            temp_nodes[v2.second].parent = temp_nodes.size(); // parent is new node
            size_type frq_sum = v1.first + v2.first;
            pq.push(tPII(frq_sum, temp_nodes.size()));
            temp_nodes.emplace_back(pc_node(frq_sum, 0, pc_node::undef,
                                            v1.second, v2.second));
        }
    }
};

struct huff_shape {
    template<class t_wt>
    using type = _huff_shape<t_wt>;
};


}// end namespace sdsl
#endif
