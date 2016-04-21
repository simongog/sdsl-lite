/* sdsl - succinct data structures library
    Copyright (C) 2009-2013 Simon Gog

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
/*! \file wt_blcd.hpp
 *  \brief wt_blcd.hpp contains a generic wavelet tree class.
 *  \author Simon Gog
*/
#ifndef INCLUDED_SDSL_WT_BLCD
#define INCLUDED_SDSL_WT_BLCD

#include "wt_pc.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

// forward declaration
struct balanced_shape;

//! A balanced wavelet tree.
/*!
 * \par Space complexity
 *      \f$\Order{n\log|\Sigma| + 2|\Sigma|\log n}\f$ bits, where \f$n\f$ is
 *      the size of the vector the wavelet tree was build for.
 *
 * \tparam t_bitvector   Underlying bitvector structure.
 * \tparam t_rank        Type of the support structure for rank on pattern `1`.
 * \tparam t_select      Type of the support structure for select on pattern `1`.
 * \tparam t_select_zero Type of the support structure for select on pattern `0`.
 *
 * \par Reference
 *  Roberto Grossi, Ankur Gupta, Jeffrey Scott Vitter:
 *  High-order entropy-compressed text indexes.
 *  Proceedings of the 14th Annual ACM-SIAM Symposium on
 *  Discrete Algorithms (SODA 2003).
 *
 * @ingroup wt
 */
template<class t_bitvector   = bit_vector,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select_one  = typename t_bitvector::select_1_type,
         class t_select_zero = typename t_bitvector::select_0_type,
         class t_tree_strat  = byte_tree<>
         >
using wt_blcd = wt_pc<balanced_shape,
      t_bitvector,
      t_rank,
      t_select_one,
      t_select_zero,
      t_tree_strat>;

template<class t_wt>
struct _balanced_shape {
    typedef typename t_wt::size_type       size_type;
    typedef std::pair<uint64_t, uint64_t>  tPII;    // (freq, nodenr)-pair
    enum { lex_ordered = 1 };

    template<class t_rac>
    static void
    construct_tree(t_rac& C, std::vector<pc_node>& temp_nodes) {
        size_type c = 0;
        std::vector<uint64_t> symbols;
        std::for_each(std::begin(C), std::end(C), [&](decltype(*std::begin(C)) &freq) {
            if (freq > 0) {
                symbols.push_back(c);
            }
            ++c;
        });
        uint64_t sigma = symbols.size();
        if (sigma > 0) {
            _construct_tree(pc_node::undef, symbols, 0, sigma, C, temp_nodes);
            pc_node root = temp_nodes[0];
            for (uint64_t i=1; i < temp_nodes.size(); ++i) {
                temp_nodes[i-1] = temp_nodes[i];
                temp_nodes[i-1].parent = (temp_nodes[i-1].parent+temp_nodes.size()-1)%temp_nodes.size();
                temp_nodes[i-1].child[0] -= (temp_nodes[i-1].child[0] != pc_node::undef);
                temp_nodes[i-1].child[1] -= (temp_nodes[i-1].child[1] != pc_node::undef);
            }
            root.child[0] -= (root.child[0] != pc_node::undef);
            root.child[1] -= (root.child[1] != pc_node::undef);
            temp_nodes[temp_nodes.size()-1] = root;
        }
    }

    // recursive construct_tree method returns node frequency and node pointer
    template<class t_rac>
    static tPII
    _construct_tree(uint64_t parent,
                    const std::vector<uint64_t>& symbols,
                    uint64_t lb,
                    uint64_t sigma,
                    const t_rac& C, std::vector<pc_node>& temp_nodes) {
        if (sigma == 1) {
            uint64_t freq = C[symbols[lb]];
            temp_nodes.emplace_back(pc_node(freq, symbols[lb], parent, pc_node::undef, pc_node::undef));
            return tPII(freq, temp_nodes.size()-1);
        } else {
            temp_nodes.emplace_back(pc_node(0, 0, parent, pc_node::undef, pc_node::undef));
            uint64_t node_id = temp_nodes.size()-1;
            uint64_t l_sigma = (sigma+1)/2;
            tPII freq_nptr_0 = _construct_tree(node_id, symbols, lb, l_sigma, C, temp_nodes);
            tPII freq_nptr_1 = _construct_tree(node_id, symbols, lb+l_sigma, sigma-l_sigma, C, temp_nodes);
            uint64_t freq = freq_nptr_0.first + freq_nptr_1.first;
            temp_nodes[node_id].freq = freq;
            temp_nodes[node_id].child[0] = freq_nptr_0.second;
            temp_nodes[node_id].child[1] = freq_nptr_1.second;
            return tPII(freq, node_id);
        }
    }

};

struct balanced_shape {
    template<class t_wt>
    using type = _balanced_shape<t_wt>;
};

}// end namespace sdsl
#endif
