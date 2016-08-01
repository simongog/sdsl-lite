/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

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
/*! \file hybrid_k2_tree.hpp
    \brief hybrid_k2_tree.hpp contains a compact hybrid k^2-tree.
    \author Jan Broß, based on the k2 treap code of Simon Gog
*/
#ifndef INCLUDED_SDSL_HYBRID_COMPRESSED_K2_TREE
#define INCLUDED_SDSL_HYBRID_COMPRESSED_K2_TREE

#include "vectors.hpp"
#include "bits.hpp"
#include "k2_tree_helper.hpp"
#include "k2_tree_base.hpp"
#include "../../external/dacs/include/dacs.h"
#include "k2_tree_vocabulary.h"
#include <tuple>
#include <algorithm>
#include <stxxl/vector>
#include <climits>
#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <gtest/gtest_prod.h>

//! Namespace for the succinct data structure library.
namespace sdsl {

/*! A hybrid k2 tree implementation
 *  \par References
 *       [1] Nieves R. Brisaboa, Susana Ladra, and Gonzalo Navarro:
 *           Compact representation of Web graphs with extended functionality.
 *           Inf. Syst. 39 (January 2014), 152-174. DOI=http://dx.doi.org/10.1016/j.is.2013.08.003
 */
    template<uint8_t t_k_l_1,
            uint8_t t_k_l_1_size,
            uint8_t t_k_l_2,
            uint8_t t_k_leaves,
            typename t_lev=bit_vector,
            typename t_leaf=bit_vector,
            typename t_rank=typename t_lev::rank_1_type>
    class k2_tree_hybrid_compressed : public k2_tree_base<t_lev,t_leaf,t_rank> {
        static_assert(t_k_l_1 > 1, "t_k has to be larger than 1.");
        static_assert(t_k_l_1 <= 16, "t_k has to be smaller than 17.");
        static_assert(t_k_l_2 > 1, "t_k has to be larger than 1.");
        static_assert(t_k_l_2 <= 16, "t_k has to be smaller than 17.");
        static_assert(t_k_leaves > 1, "t_k has to be larger than 1.");
        static_assert(t_k_leaves <= 16, "t_k has to be smaller than 17.");

    public:

        std::shared_ptr<Vocabulary> m_vocabulary;
        sFTRep *m_comp_leaves;

        k2_tree_hybrid_compressed(std::vector<t_lev> levels, std::vector<t_rank> level_rank, sFTRep *leafs,
                                  std::shared_ptr<Vocabulary> vocabulary, uint8_t tree_height, uint64_t max_element,
                                  size_type size, uint8_t access_shortcut_size, bit_vector access_shortcut,
                                  rank_support_v<01, 2> access_shortcut_rank_01,
                                  int_vector<1>::select_1_type access_shortcut_select_1):  m_vocabulary(vocabulary), m_comp_leaves(leafs) {

            this->m_levels.resize(levels.size());
            for (uint64_t j = 0; j < this->m_levels.size(); ++j) {
                this->m_levels[j].swap(levels[j]);
            }

            this->m_levels_rank.resize(this->m_levels.size());
            for (uint64_t j = 0; j < this->m_levels.size(); ++j) {
                util::swap_support(this->m_levels_rank[j], level_rank[j], &this->m_levels[j], &(levels[j]));
            }

            this->m_tree_height = tree_height;
            this->m_max_element = max_element;
            this->m_size = size;
            this->m_access_shortcut_size = access_shortcut_size;

            this->m_access_shortcut.swap(access_shortcut);
            util::swap_support(this->m_access_shortcut_rank_01_support, access_shortcut_rank_01,
                               &this->m_access_shortcut, &access_shortcut);
            util::swap_support(this->m_access_shortcut_select_1_support, access_shortcut_select_1,
                               &this->m_access_shortcut, &access_shortcut);

        }

        typedef stxxl::VECTOR_GENERATOR<std::pair<uint32_t, uint32_t>>::result stxxl_32bit_pair_vector;
        typedef stxxl::VECTOR_GENERATOR<std::pair<uint64_t, uint64_t>>::result stxxl_64bit_pair_vector;
        using k2_tree_base<t_lev,t_leaf,t_rank>::operator=;
        using k2_tree_base<t_lev,t_leaf,t_rank>::operator==;

        k2_tree_hybrid_compressed() = default;

        /*k2_tree_hybrid_compressed(const k2_tree_hybrid_compressed &tr) {
            *this = tr;
        }

        k2_tree_hybrid_compressed(k2_tree_hybrid_compressed &&tr) {
            *this = std::move(tr);
        }*/


        inline uint8_t get_k(uint8_t level) const {
            if (level < t_k_l_1_size) {
                return t_k_l_1;
            } else if (level < this->m_tree_height - 1) {
                return t_k_l_2;
            } else {
                return t_k_leaves;
            }
        }

        /**
         * Checks whether link from p = link.first to q = link.second is present i.e. matrix entry a_pq = 1
         */
        template<typename t_x, typename t_y>
        bool check_link(std::pair<t_x,t_y> link) const{

            //Patological case happening e.g. when using k2part
            if (this->m_tree_height == 0){
                return false;
            }

            return check_link_internal(0, this->m_max_element, link.first, link.second, 0);
        }

        virtual size_type serialize(std::ostream &out, structure_tree_node *v, std::string name) const override {
            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = k2_tree_base<t_lev,t_leaf,t_rank>::serialize(out, child, name);
            written_bytes += m_vocabulary->serialize(out,child,name);
            SaveFT(out, m_comp_leaves);
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        virtual void load(std::istream &in) override {
            k2_tree_base<t_lev,t_leaf,t_rank>::load(in);
            m_vocabulary->load(in);
            m_comp_leaves = LoadFT(in);
        }

    private:
        /**
        * Returns word containing the bit at the given position
        * It access the corresponding word in the DAC.
        *
        * @param pos Position in the complete sequence of bit of the last level.
        * @return Pointer to the first position of the word.
        */
        inline bool is_leaf_bit_set(uint64_t pos, uint8_t leafK) const {
            uint iword = accessFT(m_comp_leaves, pos/(leafK*leafK));
            const uchar * word = m_vocabulary->get(iword);
            bool bitSet = ((word[pos/kUcharBits] >> (pos%kUcharBits)) & 1);
            return bitSet;
        }

        /**
        * Returns word containing the bit at the given position
        * It access the corresponding word in the DAC.
        *
        * @param pos Position in the complete sequence of bit of the last level.
        * @return Pointer to the first position of the word.
        */
        template<typename t_x>
        inline void  check_leaf_bits_direct(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> & result) const {
            uint iword = accessFT(m_comp_leaves, pos/(leafK*leafK));
            const uchar * word = m_vocabulary->get(iword);
            for (int i = 0; i < leafK; ++i) {
                if ((word[(pos+i)/kUcharBits] >> ((pos+i)%kUcharBits)) & 1){
                    result.push_back(i+result_offset);
                }
            }
        }

        /**
        * Returns word containing the bit at the given position
        * It access the corresponding word in the DAC.
        *
        * @param pos Position in the complete sequence of bit of the last level.
        * @return Pointer to the first position of the word.
        */
        template<typename t_x>
        inline void  check_leaf_bits_inverse(int64_t pos, t_x result_offset, uint8_t leafK, std::vector<t_x> & result) const {
            uint iword = accessFT(m_comp_leaves, pos/(leafK*leafK));
            const uchar * word = m_vocabulary->get(iword);
            for (int i = 0; i < get_k(leafK); ++i) {
                if ((word[(pos+i*get_k(leafK))/kUcharBits] >> ((pos+i*leafK)%kUcharBits)) & 1){
                    result.push_back(i+result_offset);
                }
            }
        }
        /**
         * Checks wether the edge p-->q exists recursively
         * @param level
         *  current level, initialy 0
         * @param p
         *  source_node
         * @param q
         *  target_node
         * @param index
         *  contains the index of the first child of the previous node, initially set to 0
         * @return
         */
        template<typename t_x, typename t_y>
        bool check_link_internal(uint8_t level, uint64_t n, t_x p, t_y q, int64_t index) const {
            using namespace k2_treap_ns;

            const uint8_t k = get_k(level);

            uint64_t current_submatrix_size = n/k;
            int64_t y = index + k * (p/current_submatrix_size) + (q/current_submatrix_size);

            if (this->is_leaf_level(level)){
                return is_leaf_bit_set(y,get_k(this->m_tree_height-1));
            } else if (this->m_levels[level][y]) {
                return check_link_internal(level+1, current_submatrix_size, p % current_submatrix_size, q % current_submatrix_size, this->get_child_index(0, y, level));
            } else {
                return false;
            }
        }

        /**
        * Variant of direct_links2_internal using a queue
        * @param source_id
        * @param result
        */
        template<typename t_x>
        void  direct_links2_internal_queue(t_x source_id, std::vector<t_x> &result) const override{
            using namespace k2_treap_ns;
            //n, level, source_id, column_offset, index
            std::queue<std::tuple<uint64_t, uint8_t, t_x,t_x,int64_t>> queue;
            queue.push(std::make_tuple(this->m_max_element, 0, source_id, (t_x) 0, 0));

            while (!queue.empty()){
                auto current_element = queue.front();
                uint64_t n = std::get<0>(current_element);
                uint8_t level = std::get<1>(current_element);
                t_x source_id = std::get<2>(current_element);
                t_x column_offset = std::get<3>(current_element);
                int64_t index = std::get<4>(current_element);
                queue.pop();

                const uint8_t k = get_k(level);

                uint64_t submatrix_size = n/k;
                int64_t y = index*k*k + k *(source_id/submatrix_size);

                if (this->is_leaf_level(level)){
                    check_leaf_bits_direct(y,column_offset,get_k(level),result);
                } else { //internal node
                    for (uint j = 0; j < k; ++j) {
                        if (this->m_levels[level][y+j] == 1) {
                            queue.push(std::make_tuple(submatrix_size, level + 1, t_x(source_id % submatrix_size),
                                                       t_x(column_offset + submatrix_size * j), this->m_levels_rank[level](y+j)));
                        }
                    }
                }
            }
        }

        /**
        * Variant of inverse_links2_internal using a queue
        * @param source_id
        * @param result
        */
        template<typename t_x>
        void inverse_links2_internal_queue(t_x source_id, std::vector<t_x> &result) const {
            using namespace k2_treap_ns;
            //n, level, source_id, column_offset, index
            std::queue<std::tuple<uint64_t,uint8_t,t_x,t_x,int64_t>> queue;
            queue.push(std::make_tuple(this->m_max_element, 0, source_id, (t_x) 0, 0));

            while (!queue.empty()) {
                auto current_element = queue.front();
                t_x n = std::get<0>(current_element);
                uint8_t level = std::get<1>(current_element);
                t_x source_id = std::get<2>(current_element);
                t_x row_offset = std::get<3>(current_element);
                int64_t index = std::get<4>(current_element);
                queue.pop();

                const uint8_t k = get_k(level);

                uint64_t submatrix_size = n / k;
                int64_t y = index * k * k + (source_id / submatrix_size);

                if (this->is_leaf_level(level)) {
                    check_leaf_bits_inverse(y, row_offset,get_k(level),result);
                } else { //internal node
                    for (int j = 0; j < k; ++j) {
                        if (this->m_levels[level][y + (j * k)]){
                            queue.push(std::make_tuple(submatrix_size, level + 1, t_x(source_id % submatrix_size),
                                                       t_x(row_offset + submatrix_size * j), this->m_levels_rank[level](y + (j * k))));
                        }
                    }
                }
            }
        }

    };
}
#endif

