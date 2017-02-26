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
/*! \file bp_support_sada.hpp
    \brief bp_support_sada.hpp contains an implementation of a balanced
     parentheses support structure proposed by Kunihiko Sadakane.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_BP_SUPPORT_SADA
#define INCLUDED_SDSL_BP_SUPPORT_SADA

#include "int_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "bp_support_algorithm.hpp"
#include "fast_cache.hpp"
#include <stack>
#include <map>
#include <set>
#include <utility>
#include <stdexcept>
#ifndef NDEBUG
#include <algorithm>
#endif
#include <iostream>

namespace sdsl
{

//! A class that provides support for bit_vectors that represent a BP sequence.
/*! This data structure supports the following operations:
 *   - find_open
 *   - find_close
 *   - enclose
 *   - double_enclose
 *   - rank
 *   - select
 *   - excess
 *   - rr_enclose
 *  An opening parenthesis in the balanced parentheses sequence is represented by a 1 in the bit_vector
 *  and a closing parenthesis by a 0.
 *
 *  \tparam t_sml_blk The size of the small blocks. Denoted as `s` in Sadakane's paper.
 *  \tparam t_med_deg Number of small blocks that a medium block contains. Denoted as `l` in Sadakane's paper.
 *  \tparam t_rank    Type of rank support used for the underlying bitvector.
 *  \tparam t_select  Type of select support used for the underlying bitvector.
 *
 *  \par References
 *      - Kunihiko Sadakane:
 *        The Ultimate Balanced Parentheses
 *        Technical Report 2008.
 *      - Kunihiko Sadakane, Gonzalo Navarro:
 *        Fully-Functional Succinct Trees.
 *        SODA 2010: 134-149
 *
 *  @ingroup bps
 */
template<uint32_t t_sml_blk = 256,
         uint32_t t_med_deg = 32,
         class t_rank       = rank_support_v5<>,
         class t_select     = select_support_mcl<> >
class bp_support_sada
{
    public:
        typedef bit_vector::size_type       size_type;
        typedef bit_vector::difference_type difference_type;
        typedef int_vector<>                sml_block_array_type;
        typedef int_vector<>                med_block_array_type;
        typedef t_rank                      rank_type;
        typedef t_select                    select_type;
    private:
        static_assert(0 < t_sml_blk, "bp_support_sada: t_sml_blk should be greater than 0!");
        const bit_vector* m_bp        = nullptr;   // the supported balanced parentheses sequence as bit_vector
        rank_type         m_bp_rank;   // RS for the BP sequence => see excess() and rank()
        select_type       m_bp_select; // SS for the BP sequence => see select()

        sml_block_array_type  m_sml_block_min_max;
        med_block_array_type  m_med_block_min_max;

        size_type m_size             = 0; // number of supported parentheses
        size_type m_sml_blocks       = 0; // number of small sized blocks
        size_type m_med_blocks       = 0; // number of medium sized blocks
        size_type m_med_inner_blocks = 0; // number of inner nodes in the min max tree of the medium sized blocks
//#define USE_CACHE
#ifdef USE_CACHE
        mutable fast_cache find_close_cache;
        mutable fast_cache find_open_cache;
        mutable fast_cache select_cache;
#endif

        void copy(const bp_support_sada& bp_support)
        {
            m_bp        = bp_support.m_bp;
            m_bp_rank   = bp_support.m_bp_rank;
            m_bp_rank.set_vector(m_bp);
            m_bp_select = bp_support.m_bp_select;
            m_bp_select.set_vector(m_bp);

            m_sml_block_min_max = bp_support.m_sml_block_min_max;
            m_med_block_min_max = bp_support.m_med_block_min_max;

            m_size             = bp_support.m_size;
            m_sml_blocks       = bp_support.m_sml_blocks;
            m_med_blocks       = bp_support.m_med_blocks;
            m_med_inner_blocks = bp_support.m_med_inner_blocks;
        }

        inline static size_type sml_block_idx(size_type i)
        {
            return i/t_sml_blk;
        }

        inline static size_type med_block_idx(size_type i)
        {
            return i/(t_sml_blk*t_med_deg);
        }

        inline static bool is_root(size_type v)
        {
            return v==0;
        }

        inline static bool is_left_child(size_type v)
        {
            assert(!is_root(v));
            return v%2;
        }

        inline static bool is_right_child(size_type v)
        {
            assert(!is_root(v));
            return !(v%2);
        }

        inline static size_type parent(size_type v)
        {
            assert(!is_root(v));
            return (v-1)/2;
        }

        inline static size_type left_child(size_type v)
        {
            return 2*v+1;
        }

        inline static size_type right_child(size_type v)
        {
            return 2*v+2;
        }

        inline bool node_exists(size_type v)const
        {
            return v < (m_med_inner_blocks + m_med_blocks);
        }

        inline static size_type right_sibling(size_type v)
        {
            return ++v;
        }

        inline static size_type left_sibling(size_type v)
        {
            return --v;
        }

        inline bool is_leaf(size_type v)const
        {
            return v >= m_med_inner_blocks;
        }

        inline difference_type min_value(size_type v)const
        {
            return m_size-((difference_type)m_med_block_min_max[2*v]);
        }

        inline difference_type max_value(size_type v)const
        {
            return m_med_block_min_max[2*v+1]-m_size;
        }

        inline difference_type sml_min_value(size_type sml_block)const
        {
            return (1 - ((difference_type)m_sml_block_min_max[sml_block<<1]));
        }

        inline difference_type sml_max_value(size_type sml_block)const
        {
            return (difference_type)m_sml_block_min_max[(sml_block<<1)+1] - 1;
        }

        void print_node(size_type v)const
        {
            std::cout<< "v = "<< v <<"  (" << min_value(v)
                     << ", " << max_value(v) << ")" ;
            if (is_leaf(v)) {
                std::cout<<" range: ["<<(v-m_med_inner_blocks)*t_med_deg* t_sml_blk
                         << ","<<(v-m_med_inner_blocks+1)*t_med_deg* t_sml_blk-1<<"]";
            }
            std::cout<< std::endl;
        }


        //! Calculate the min parenthesis \f$j>i\f$ with \f$excess(j)=excess(i)+rel\f$
        /*! \param i   The index of a parenthesis in the supported sequence.
         *  \param rel The excess difference to the excess value of parentheses \f$i\f$.
         *  \return    If there exists a parenthesis \f$ j>i\f$ with
         *            \f$ excess(j) = excess(i)+rel \f$, \f$j\f$ is returned
         *                otherwise size().
         */
        size_type fwd_excess(size_type i, difference_type rel)const
        {
            size_type j;
            // (1) search the small block for the answer
            if ((j = near_fwd_excess(*m_bp, i+1, rel, t_sml_blk)) > i) {
                return j;
            }
            difference_type desired_excess = excess(i)+rel;
            // (2) scan the small blocks of the current median block for an answer
            if ((j = fwd_excess_in_med_block(sml_block_idx(i)+1, desired_excess)) != size()) {
                return j;
            }
            // (3) search the min-max tree of the medium blocks for the right med block
            if (med_block_idx(i) == m_med_blocks)  // if we are already in the last medium block => we are done
                return size();
            size_type v    = m_med_inner_blocks + med_block_idx(i);
            // (3 a) go up the tree
            while (!is_root(v)) {
                if (is_left_child(v)) { // if the node is a left child
                    v = right_sibling(v); // choose right sibling
                    if (min_value(v) <= desired_excess and desired_excess <= max_value(v))  // found solution
                        break;
                }
                v = parent(v); // choose parent
            }
            // (3 b) go down the tree
            if (!is_root(v)) { // found solution for the query
                while (!is_leaf(v)) {
                    v = left_child(v); // choose left child
                    if (!(min_value(v) <= desired_excess and desired_excess <= max_value(v))) {
                        v = right_sibling(v); // choose right child == right sibling of the left child
                        assert((min_value(v) <= desired_excess and desired_excess <= max_value(v)));
                    }
                }
                return fwd_excess_in_med_block((v-m_med_inner_blocks)*t_med_deg, desired_excess);
            }
            // no solution found
            return size();
        }

        //! Calculate the maximal parenthesis \f$ j<i \f$ with \f$ excess(j) = excess(i)+rel \f$
        /*! \param i     The index of a parenthesis in the supported sequence.
         *  \param rel  The excess difference to the excess value of parenthesis \f$i\f$.
         *  \return     If there exists a parenthesis \f$i<j\f$ with \f$ excess(j) = excess(i)+rel\f$, \f$j\f$ is returned
         *                otherwise size().
         */
        size_type bwd_excess(size_type i, difference_type rel)const
        {
            size_type j;
            if (i == 0) {
                return rel == 0 ? -1 : size();
            }
            // (1) search the small block for the answer
            if ((j = near_bwd_excess(*m_bp, i-1, rel, t_sml_blk)) < i or j == (size_type)-1) {
                return j;
            }
            difference_type desired_excess = excess(i)+rel;
            // (2) scan the small blocks of the current median block for an answer
            if ((j = bwd_excess_in_med_block(sml_block_idx(i)-1, desired_excess)) != size()) {
                return j;
            }
            // (3) search the min-max tree of the medium blocks for the right med block
            if (med_block_idx(i) == 0) { // if we are already in the first medium block => we are done
                if (desired_excess == 0)
                    return -1;
                return size();
            }
            size_type v    = m_med_inner_blocks + med_block_idx(i);
            // (3 a) go up the tree
            while (!is_root(v)) {
                if (is_right_child(v)) { // if the node is a right child
                    v = left_sibling(v); // choose left sibling
                    if (min_value(v) <= desired_excess and desired_excess <= max_value(v))  // found solution
                        break;
                }
                v = parent(v); // choose parent
            }
            // (3 b) go down the tree
            if (!is_root(v)) { // found solution for the query
                while (!is_leaf(v)) {
                    v = right_child(v); // choose  child
                    if (!(min_value(v) <= desired_excess and desired_excess <= max_value(v))) {
                        v = left_sibling(v); // choose left child == left sibling of the right child
                        assert((min_value(v) <= desired_excess and desired_excess <= max_value(v)));
                    }
                }
                return bwd_excess_in_med_block((v-m_med_inner_blocks)*t_med_deg+(t_med_deg-1), desired_excess);
            } else if (desired_excess == 0) {
                return -1;
            }
            // no solution found
            return size();
        }

        //! Calculate the maximal parentheses \f$ j \leq sml_block_idx\cdot t_sml_blk+(t_sml_blk-1) \f$ with \f$ excess(j)=desired\_excess \f$
        size_type bwd_excess_in_med_block(size_type sml_block_idx, difference_type desired_excess)const
        {
            // get the first small block in the medium block right to the current med block
            size_type first_sml_block_in_med_block = (med_block_idx(sml_block_idx*t_sml_blk))*t_med_deg;

            while ((sml_block_idx+1) and sml_block_idx >= first_sml_block_in_med_block) {
                difference_type ex         = (sml_block_idx == 0) ? 0 : excess(sml_block_idx*t_sml_blk-1);
                difference_type min_ex     = ex + (1 - ((difference_type)m_sml_block_min_max[2*sml_block_idx]));
                difference_type max_ex    = ex + (m_sml_block_min_max[2*sml_block_idx+1] - 1);

                if (min_ex <= desired_excess and desired_excess <= max_ex) {
                    size_type j = near_bwd_excess(*m_bp, (sml_block_idx+1)*t_sml_blk-1, desired_excess-excess((sml_block_idx+1)*t_sml_blk), t_sml_blk);
                    return j;
                }
                --sml_block_idx;
            }
            if (sml_block_idx == 0 and desired_excess == 0)
                return -1;
            return size();
        }

        //! Calculate the minimal parentheses \f$ j \geq sml_block_idx\cdot t_sml_blk \f$ with \f$ excess(j)=desired\_excess \f$
        size_type fwd_excess_in_med_block(size_type sml_block_idx, difference_type desired_excess)const
        {
            // get the first small block in the medium block right to the current med block
            size_type first_sml_block_nr_in_next_med_block = (med_block_idx(sml_block_idx*t_sml_blk)+1)*t_med_deg;
            if (first_sml_block_nr_in_next_med_block > m_sml_blocks)
                first_sml_block_nr_in_next_med_block = m_sml_blocks;

            assert(sml_block_idx > 0);
            while (sml_block_idx < first_sml_block_nr_in_next_med_block) {
                difference_type ex         = excess(sml_block_idx*t_sml_blk-1);
                difference_type min_ex     = ex + (1 - ((difference_type)m_sml_block_min_max[2*sml_block_idx]));
                difference_type max_ex    = ex + m_sml_block_min_max[2*sml_block_idx+1] - 1;
                if (min_ex <= desired_excess and desired_excess <= max_ex) {
                    size_type j = near_fwd_excess(*m_bp, sml_block_idx*t_sml_blk, desired_excess-ex, t_sml_blk);
                    return j;
                }
                ++sml_block_idx;
            }
            return size();
        }

    public:
        const rank_type&            bp_rank           = m_bp_rank;           //!< RS for the underlying BP sequence.
        const select_type&          bp_select         = m_bp_select;         //!< SS for the underlying BP sequence.
        const sml_block_array_type& sml_block_min_max = m_sml_block_min_max; //!< Small blocks array. Rel. min/max for the small blocks.
        const med_block_array_type& med_block_min_max = m_med_block_min_max; //!< Array containing the min max tree of the medium blocks.

        bp_support_sada() {}

        //! Constructor
        explicit bp_support_sada(const bit_vector* bp): m_bp(bp),
            m_size(bp==nullptr?0:bp->size()),
            m_sml_blocks((m_size+t_sml_blk-1)/t_sml_blk),
            m_med_blocks((m_size+t_sml_blk* t_med_deg-1)/(t_sml_blk* t_med_deg)),
            m_med_inner_blocks(0)
        {
            if (bp == nullptr or bp->size()==0)
                return;
            // initialize rank and select
            util::init_support(m_bp_rank, bp);
            util::init_support(m_bp_select, bp);

            m_med_inner_blocks = 1;
            // m_med_inner_blocks = (next power of 2 greater than or equal to m_med_blocks)-1
            while (m_med_inner_blocks < m_med_blocks) {
                m_med_inner_blocks <<= 1; assert(m_med_inner_blocks!=0);
            }
            --m_med_inner_blocks;
            assert((m_med_inner_blocks == 0) or(m_med_inner_blocks%2==1));

            m_sml_block_min_max = int_vector<>(2*m_sml_blocks, 0, bits::hi(t_sml_blk+2)+1);
            m_med_block_min_max = int_vector<>(2*(m_med_blocks+m_med_inner_blocks), 0, bits::hi(2*m_size+2)+1);

            // calculate min/max excess values of the small blocks and medium blocks
            difference_type min_ex = 1, max_ex = -1, curr_rel_ex = 0, curr_abs_ex = 0;
            for (size_type i=0; i < m_size; ++i) {
                if ((*bp)[i])
                    ++curr_rel_ex;
                else
                    --curr_rel_ex;
                if (curr_rel_ex > max_ex) max_ex = curr_rel_ex;
                if (curr_rel_ex < min_ex) min_ex = curr_rel_ex;
                if ((i+1)%t_sml_blk == 0 or i+1 == m_size) {
                    size_type sidx = i/t_sml_blk;
                    m_sml_block_min_max[2*sidx    ] = -(min_ex-1);
                    m_sml_block_min_max[2*sidx + 1] = max_ex+1;

                    size_type v = m_med_inner_blocks + sidx/t_med_deg;

                    if ((difference_type)(-(curr_abs_ex + min_ex)+m_size) > ((difference_type)m_med_block_min_max[2*v])) {
                        assert(curr_abs_ex+min_ex <= min_value(v));
                        m_med_block_min_max[2*v] = -(curr_abs_ex + min_ex)+m_size;
                    }

                    if ((difference_type)(curr_abs_ex + max_ex + m_size) > (difference_type)m_med_block_min_max[2*v + 1])
                        m_med_block_min_max[2*v + 1] = curr_abs_ex + max_ex + m_size;

                    curr_abs_ex += curr_rel_ex;
                    min_ex = 1; max_ex = -1; curr_rel_ex = 0;
                }
            }

            for (size_type v = m_med_block_min_max.size()/2 - 1; !is_root(v); --v) {
                size_type p = parent(v);
                if (min_value(v) < min_value(p))  // update minimum
                    m_med_block_min_max[2*p] = m_med_block_min_max[2*v];
                if (max_value(v) > max_value(p))  // update maximum
                    m_med_block_min_max[2*p+1] = m_med_block_min_max[2*v+1];
            }
        }

        //! Copy constructor
        bp_support_sada(const bp_support_sada& bp_support)
        {
            copy(bp_support);
        }

        //! Move constructor
        bp_support_sada(bp_support_sada&& bp_support)
        {
            *this = std::move(bp_support);
        }

        //! Assignment operator
        bp_support_sada& operator=(bp_support_sada&& bp_support)
        {
            if (this != &bp_support) {
                m_bp        = std::move(bp_support.m_bp);
                m_bp_rank   = std::move(bp_support.m_bp_rank);
                m_bp_rank.set_vector(m_bp);
                m_bp_select = std::move(bp_support.m_bp_select);
                m_bp_select.set_vector(m_bp);

                m_sml_block_min_max = std::move(bp_support.m_sml_block_min_max);
                m_med_block_min_max = std::move(bp_support.m_med_block_min_max);

                m_size             = std::move(bp_support.m_size);
                m_sml_blocks       = std::move(bp_support.m_sml_blocks);
                m_med_blocks       = std::move(bp_support.m_med_blocks);
                m_med_inner_blocks = std::move(bp_support.m_med_inner_blocks);
            }
            return *this;
        }

        //! Swap method
        /*! Swaps the content of the two data structure.
         *  You have to use set_vector to adjust the supported bit_vector.
         *  \param bp_support Object which is swapped.
         */
        void swap(bp_support_sada& bp_support)
        {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            m_bp_rank.swap(bp_support.m_bp_rank);
            m_bp_select.swap(bp_support.m_bp_select);

            m_sml_block_min_max.swap(bp_support.m_sml_block_min_max);
            m_med_block_min_max.swap(bp_support.m_med_block_min_max);

            std::swap(m_size, bp_support.m_size);
            std::swap(m_sml_blocks, bp_support.m_sml_blocks);
            std::swap(m_med_blocks, bp_support.m_med_blocks);
            std::swap(m_med_inner_blocks, bp_support.m_med_inner_blocks);
        }

        //! Assignment operator
        bp_support_sada& operator=(const bp_support_sada& bp_support)
        {
            if (this != &bp_support) {
                copy(bp_support);
            }
            return *this;
        }
        /*
                bp_support_sada& operator=(bp_support_sada&& bps)
                {
                    this->swap(bps);
                    return *this;
                }
        */
        void set_vector(const bit_vector* bp)
        {
            m_bp = bp;
            m_bp_rank.set_vector(bp);
            m_bp_select.set_vector(bp);
        }

        /*! Calculates the excess value at index i.
         * \param i The index of which the excess value should be calculated.
         */
        inline difference_type excess(size_type i)const
        {
            return (m_bp_rank(i+1)<<1)-i-1;
        }

        /*! Returns the number of opening parentheses up to and including index i.
         * \pre{ \f$ 0\leq i < size() \f$ }
         */
        size_type rank(size_type i)const
        {
            return m_bp_rank(i+1);
        }

        /*! Returns the index of the i-th opening parenthesis.
         * \param i Number of the parenthesis to select.
         * \pre{ \f$1\leq i < rank(size())\f$ }
         * \post{ \f$ 0\leq select(i) < size() \f$ }
         */
        size_type select(size_type i)const
        {
#ifdef USE_CACHE
            size_type a = 0;
            if (select_cache.exists(i, a)) {
                return a;
            } else {
                a = m_bp_select(i);
                select_cache.write(i, a);
                return a;
            }
#endif
            return m_bp_select(i);
        }

        /*! Calculate the index of the matching closing parenthesis to the parenthesis at index i.
         * \param i Index of an parenthesis. 0 <= i < size().
         * \return * i, if the parenthesis at index i is closing,
         *         * the position j of the matching closing parenthesis, if a matching parenthesis exists,
         *         * size() if no matching closing parenthesis exists.
         */
        size_type find_close(size_type i)const
        {
            assert(i < m_size);
            if (!(*m_bp)[i]) {// if there is a closing parenthesis at index i return i
                return i;
            }
#ifdef USE_CACHE
            size_type a = 0;
            if (find_close_cache.exists(i, a)) {
                return a;
            } else {
                a = fwd_excess(i, -1);
                find_close_cache.write(i, a);
                return a;
            }
#endif
            return fwd_excess(i, -1);
        }

        //! Calculate the matching opening parenthesis to the closing parenthesis at position i
        /*! \param i Index of a closing parenthesis.
          * \return * i, if the parenthesis at index i is closing,
          *         * the position j of the matching opening parenthesis, if a matching parenthesis exists,
          *         * size() if no matching closing parenthesis exists.
          */
        size_type find_open(size_type i)const
        {
            assert(i < m_size);
            if ((*m_bp)[i]) {// if there is a opening parenthesis at index i return i
                return i;
            }
#ifdef USE_CACHE
            size_type a = 0;
            if (find_open_cache.exists(i, a)) {
                return a;
            } else {
                size_type bwd_ex = bwd_excess(i,0);
                if (bwd_ex == size())
                    a = size();
                else
                    a = bwd_ex+1;
                find_open_cache.write(i, a);
                return a;
            }
#endif
            size_type bwd_ex = bwd_excess(i,0);
            if (bwd_ex == size())
                return size();
            else
                return bwd_ex+1;
        }

        //! Calculate the index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i.
        /*! \param i Index of an opening parenthesis.
         *  \return The index of the opening parenthesis corresponding to the closest matching parenthesis pair enclosing i,
         *          or size() if no such pair exists.
         */
        size_type enclose(size_type i)const
        {
            assert(i < m_size);
            if (!(*m_bp)[i]) { // if there is closing parenthesis at position i
                return find_open(i);
            }
            size_type bwd_ex = bwd_excess(i, -2);
            if (bwd_ex == size())
                return size();
            else
                return bwd_ex+1;
        }

        //! The range restricted enclose operation for parentheses pairs \f$(i,\mu(i))\f$ and \f$(j,\mu(j))\f$.
        /*! \param i First opening parenthesis.
         *  \param j Second opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The smallest index, say k, of an opening parenthesis such that findclose(i) < k < j and
         *  findclose(j) < findclose(k). If such a k does not exists, restricted_enclose(i,j) returns size().
         * \par Time complexity
         *     \f$ \Order{block\_size} \f$
         */
        size_type rr_enclose(const size_type i, const size_type j)const
        {
            assert(j < m_size);
            assert((*m_bp)[i]==1 and(*m_bp)[j]==1);
            const size_type mip1 = find_close(i)+1;
            if (mip1 >= j)
                return size();
            return rmq_open(mip1, j);
        }

        /*! Search the interval [l,r-1] for an opening parenthesis, say i, such that find_close(i) >= r.
         * \param l The left end (inclusive) of the interval to search for the result.
         * \param r The right end (exclusive) of the interval to search for the result.
         * \return The minimal opening parenthesis i with \f$ \ell \leq i < r \f$ and \f$ find_close(i) \geq r \f$;
         *         if no such i exists size() is returned.
         * \par Time complexity
         *     \f$ \Order{block\_size} \f$
         */
        size_type rmq_open(const size_type l, const size_type r)const
        {
            assert(r < m_bp->size());
            if (l >= r)
                return size();
            size_type res = rmq(l, r-1);
            assert(res>=l and res<=r-1);
            if ((*m_bp)[res] == 1) { // The parenthesis with minimal excess is opening
                assert(find_close(res) >= r);
                return res;
            } else {
                res = res+1; // go to the next parenthesis to the right
                if (res < r) { // The parenthesis with minimal excess if closing and the next opening parenthesis is less than r
                    assert((*m_bp)[res] == 1);
                    size_type ec = enclose(res);
                    if (ec < l or ec == size()) {
                        assert(find_close(res)>=r);
                        return res;
                    } else {
                        assert(find_close(ec)>=r);
                        return ec;
                    }
                } else if (res == r) {
                    size_type ec = enclose(res); // if m_bp[res]==0 => find_open(res), if m_bp[res]==1 => enclose(res)
                    if (ec >= l) {
                        assert(ec == size() or excess(ec)==excess(res-1));
                        return ec;
                    }
                }
            }
            return size();
        }

        size_type median_block_rmq(size_type l_sblock, size_type r_sblock, bit_vector::difference_type& min_ex)const
        {
            assert(l_sblock <= r_sblock+1);
            size_type pos_min_block = (size_type)-1;
            difference_type e = 0;
            if (l_sblock == 0) {
                if (sml_min_value(0) <= min_ex) {
                    pos_min_block = 0;
                    min_ex = sml_min_value(0);
                }
                l_sblock = 1;
            }
            for (size_type i=l_sblock; i <= r_sblock; ++i) {
                if ((e = (excess(i*t_sml_blk-1)  + sml_min_value(i))) <= min_ex) {
                    pos_min_block = i;
                    min_ex = e;
                }
            }
            return pos_min_block;
        }


        //! The range minimum query (rmq) returns the index of the parenthesis with minimal excess in the range \f$[l..r]\f$
        /*! \param l The left border of the interval \f$[l..r]\f$ (\f$l\leq r\f$).
         *  \param r The right border of the interval \f$[l..r]\f$ (\f$l \leq r\f$).
         */
        size_type rmq(size_type l, size_type r)const
        {
            assert(l<=r);
            size_type sbl = sml_block_idx(l);
            size_type sbr = sml_block_idx(r);
            difference_type min_rel_ex = 0;
            if (sbl == sbr) { // if l and r are in the same small block
                return near_rmq(*m_bp, l, r, min_rel_ex);
            } else {
                difference_type min_ex  = 0;          // current minimal excess value
                size_type        min_pos = 0;          // current min pos
                enum min_pos_type {POS, SMALL_BLOCK_POS, MEDIUM_BLOCK_POS};
                enum min_pos_type pos_type = POS;     // current
                min_pos = near_rmq(*m_bp, l, (sbl+1)*t_sml_blk-1, min_rel_ex); // scan the leftmost small block of l
                assert(min_pos >= l);
                min_ex = excess(l) + min_rel_ex;

                size_type mbl = med_block_idx(l);
                size_type mbr = med_block_idx(r);    assert(mbl <= mbr);

                size_type temp = median_block_rmq(sbl+1, std::min((mbl+1)*t_med_deg-1, sbr-1), min_ex); // scan the medium block of l
                if (temp != (size_type)-1) {
                    assert(temp* t_sml_blk >= l and temp* t_sml_blk <= r);
                    min_pos     = temp;
                    assert(min_pos >= 0  and min_pos < m_sml_blocks);
                    pos_type     = SMALL_BLOCK_POS;
                }
#if 0
                // sequential scan the medium blocks
                for (size_type v=mbl+1+m_med_inner_blocks; v < mbr + m_med_inner_blocks; ++v) {
                    assert(is_leaf(v));
                    if (min_value(v) <= min_ex) {
                        min_ex = min_value(v);
                        min_pos = v;
                        assert(min_pos-m_med_inner_blocks >= 0 and min_pos < m_med_blocks-m_med_inner_blocks);
                        pos_type = MEDIUM_BLOCK_POS;
                    }
                }
#else
                // tree search in the min max tree of the medium blocks
                if (mbr-mbl > 1) {
                    size_type v = mbl + 1 + m_med_inner_blocks;
                    size_type rcb = mbl + 1; // rightmost covered block
                    size_type h = 0; // subtree height
                    while (rcb < mbr-1) {  // go up the tree until the rightmost covered block >= mbr-1
                        if (min_value(v) <= min_ex) {
                            min_ex = min_value(v); min_pos = v; pos_type = MEDIUM_BLOCK_POS;
                        }
                        if (is_right_child(v)) { // v is a right child
                            h += 1;
                            rcb += (1ULL<<h);
                            v = right_sibling(parent(v));
                        } else { //  it is a left child
                            rcb += (1ULL<<h);
                            h += 1;
                            v = parent(v);
                        }
                    }
                    if (rcb <= mbr-1 and min_value(v) <= min_ex) {
                        min_ex = min_value(v); min_pos = v; pos_type = MEDIUM_BLOCK_POS;
                    }
                    assert(node_exists(v));
                    assert(rcb >= mbr-1);

                    while (rcb != mbr-1) { // go down the tree until the rightmost covered block = mbr-1
                        assert(h != (size_type)-1);
                        if (rcb > mbr-1) {
                            h     = h-1;
                            rcb = rcb - (1ULL<<h);
                            v = left_child(v);
                        } else { // rcb < mbr-1
                            h = h-1;
                            rcb = rcb + (1ULL<<h);
                            v = right_sibling(right_child(v));
                        }
                        if (rcb <= mbr-1 and min_value(v) <= min_ex) {
                            min_ex = min_value(v); min_pos = v; pos_type = MEDIUM_BLOCK_POS;
                        }
                    }
                    if (pos_type == MEDIUM_BLOCK_POS) {
                        while (!is_leaf(min_pos)) {
                            min_pos = right_child(min_pos);
                            if (!node_exists(min_pos) or min_value(min_pos) > min_ex)
                                min_pos = left_sibling(min_pos);
                        }
                    }
                }
#endif

                // search in the medium block of r
                temp = median_block_rmq(std::max(mbr*t_med_deg, sbl+1), sbr-1, min_ex);  // scan the medium block of r
                if (temp != (size_type)-1) {
                    assert(temp* t_sml_blk >= l and temp* t_sml_blk <= r);
                    min_pos     = temp;
                    pos_type     = SMALL_BLOCK_POS;
                }
                // search in the small block of r
                temp = near_rmq(*m_bp, sbr*t_sml_blk, r, min_rel_ex); // scan the small block of r
                if ((excess(sbr*t_sml_blk) + min_rel_ex) <= min_ex) {             // if it contains the minimum return its position
                    assert(temp>=l and temp<=r);
                    return temp;
                }
                // if the found minimum lies in a medium block => find its small block
                if (pos_type == MEDIUM_BLOCK_POS) {
                    min_pos = min_pos - m_med_inner_blocks;
                    temp = median_block_rmq(min_pos*t_med_deg, (min_pos+1)*t_med_deg-1, min_ex);
                    assert(temp != (size_type)-1);   // assert that we find a solution
                    assert(temp* t_sml_blk >= l and temp* t_sml_blk <= r);
                    min_pos = temp;
                    pos_type = SMALL_BLOCK_POS;
                }
                if (pos_type == SMALL_BLOCK_POS) {
                    min_pos = near_rmq(*m_bp, min_pos*t_sml_blk, (min_pos+1)*t_sml_blk-1, min_rel_ex);
                    assert(min_pos >=l and min_pos <= r);
                }
                return min_pos;
            }
        }

        //! The range restricted enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The minimal opening parenthesis, say k, such that \f$ findclose(i) < k < j\f$ and
         *  findclose(j) < findclose(k). If such a k does not exists, restricted_enclose(i,j) returns size().
         *  \par Time complexity
         *        \f$ \Order{size()}\f$ in the worst case.
        */
        size_type rr_enclose_naive(size_type i, size_type j)const
        {
            assert(j > i and j < m_size);
            assert((*m_bp)[i]==1 and(*m_bp)[j]==1);
            size_type mi = find_close(i); // matching parenthesis to i
            assert(mi > i and mi < j);
            assert(find_close(j) > j);
            size_type k = enclose(j);
            if (k == m_size or k < i) // there exists no opening parenthesis at position mi<k<j.
                return m_size;
            size_type kk;
            do {
                kk = k;
                k = enclose(k);
            } while (k != m_size and k > mi);
            return kk;
        }

        //! The double enclose operation
        /*! \param i Index of an opening parenthesis.
         *  \param j Index of an opening parenthesis \f$ i<j \wedge findclose(i) < j \f$.
         *  \return The maximal opening parenthesis, say k, such that \f$ k<j \wedge k>findclose(j) \f$.
         *          If such a k does not exists, double_enclose(i,j) returns size().
         */
        size_type double_enclose(size_type i, size_type j)const
        {
            assert(j > i);
            assert((*m_bp)[i]==1 and(*m_bp)[j]==1);
            size_type k = rr_enclose(i, j);
            if (k == size())
                return enclose(j);
            else
                return enclose(k);
        }

        //! Return the number of zeros which proceed position i in the balanced parentheses sequence.
        /*! \param i Index of an parenthesis.
         */
        size_type preceding_closing_parentheses(size_type i)const
        {
            assert(i < m_size);
            if (!i) return 0;
            size_type ones = m_bp_rank(i);
            if (ones) { // ones > 0
                assert(m_bp_select(ones) < i);
                return i - m_bp_select(ones) - 1;
            } else {
                return i;
            }
        }

		//! Returns the level ancestor of the node i.
		/*! \param i The index of a parenthesis (i.e., a node).
		 *  \param d The level, i.e., which node to select on the path from the node i up to the root. 
		 *           The level d = 0 will return the node itself, d = 1 will return its parent, and so on.
		 */
		size_type level_anc(size_type i, size_type d)const
		{
			assert(i < m_size);
			size_type bwd_ex = bwd_excess(i,-d-1);
			if (bwd_ex == size())
				return size();
			else
				return bwd_ex+1;
		}

		/*! The size of the supported balanced parentheses sequence.
         * \return the size of the supported balanced parentheses sequence.
         */
        size_type size() const
        {
            return m_size;
        }

        //! Serializes the bp_support_sada to a stream.
        /*!
         * \param out The outstream to which the data structure is written.
         * \return The number of bytes written to out.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += write_member(m_sml_blocks, out, child, "sml_block_cnt");
            written_bytes += write_member(m_med_blocks, out, child, "med_block_cnt");
            written_bytes += write_member(m_med_inner_blocks, out, child, "med_inner_blocks");

            written_bytes += m_bp_rank.serialize(out, child, "bp_rank");
            written_bytes += m_bp_select.serialize(out, child, "bp_select");

            written_bytes += m_sml_block_min_max.serialize(out, child, "sml_blocks");
            written_bytes += m_med_block_min_max.serialize(out, child, "med_blocks");

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load the bp_support_sada for a bit_vector v.
        /*!
         * \param in The instream from which the data structure is read.
         * \param bp Bit vector representing a balanced parentheses sequence that is supported by this data structure.
         */
        void load(std::istream& in, const bit_vector* bp)
        {
            m_bp = bp;
            read_member(m_size, in);
            assert(m_size == bp->size());
            read_member(m_sml_blocks, in);
            read_member(m_med_blocks, in);
            read_member(m_med_inner_blocks, in);

            m_bp_rank.load(in, m_bp);
            m_bp_select.load(in, m_bp);

            m_sml_block_min_max.load(in);
            m_med_block_min_max.load(in);
        }
};

}// end namespace




#endif
