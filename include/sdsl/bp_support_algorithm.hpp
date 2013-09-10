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
/*! \file bp_support_algorithm.hpp
    \brief bp_support_algorithm.hpp contains algorithms for balanced parentheses sequences.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_BP_SUPPORT_ALGORITHM
#define INCLUDED_SDSL_BP_SUPPORT_ALGORITHM

#include "int_vector.hpp" // for bit_vector
#include <stack> // for calculate_pioneers_bitmap method
#include <map>   // for calculate_pioneers_bitmap method
#include "sorted_stack_support.hpp"


namespace sdsl
{

// This structure contains lookup tables
struct excess {
    static struct impl {
        // Given an excess value x in [-8,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // near_fwd_pos[(x+8)<<8 | w] contains the minimal position
        // p in [0..7] where the excess value x is reached, or 8
        // if x is not reached in w.
        uint8_t near_fwd_pos[(8-(-8))*256];

        // Given an excess value of x in [-8,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // near_bwd_pos[(x+8)<<8 | w] contains the maximal position
        // p in [0..7] where the excess value x is reached, or 8
        // if x is not reached in w.
        uint8_t near_bwd_pos[(8-(-8))*256];

        // Given a 8-bit word w. word_sum[w] contains the
        // excess value of w.
        int8_t word_sum[256];

        // Given a 8-bit word w. min[w] contains the
        // minimal excess value in w.
        int8_t min[256];

        // Given a 8-bit word w. min_pos_max[w] contains
        // the maximal position p in w, where min[w] is
        // reached
        int8_t min_pos_max[256];

        // Given an excess value x in [1,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // min_match_pos_packed[w]:[(x-1)*4,x*4] contains
        // the minimal position, where excess value
        // -x is reached and 9, if there is no such position.
        uint32_t min_match_pos_packed[256];

        // Given an excess value x in [1,8] and a 8-bit
        // word w interpreted as parentheses sequence.
        // max_match_pos_packed[w]:[(x-1)*4,x*4] contains
        // the maximal position, where excess value
        // -x is reached and 9, if there is no such position.
        uint32_t max_match_pos_packed[256];

        // Given a 8-bit word w. x=min_and_info[w] contains
        // the following information.
        // * [0..7] the minimum excess value in w + 8 of an opening parenthesis
        // * [8..11] the maximal position of the minimal excess value
        // * [12..15] the number of ones in the word
        // if w != 0, and 17 for w=0.
        uint16_t min_open_excess_info[256];

        impl() {
            for (int32_t x = -8; x < 8; ++x) {
                for (uint16_t w=0; w < 256; ++w) {
                    uint16_t i = (x+8)<<8|w;
                    near_fwd_pos[i] = 8;
                    int8_t p=0;
                    int8_t excess = 0;
                    do {
                        excess += 1-2*((w&(1<<p))==0);
                        if (excess == x) {
                            near_fwd_pos[i] = p;
                            break;
                        }
                        ++p;
                    } while (p < 8);

                    near_bwd_pos[i] = 8;
                    p = 7;
                    excess = 0;
                    do {
                        excess += 1-2*((w&(1<<p))>0);
                        if (excess == x) {
                            near_bwd_pos[i] = p;
                            break;
                        }
                        --p;
                    } while (p > -1);
                }
            }
            int_vector<> packed_mins(1, 0, 32);
            int_vector<> packed_maxs(1, 0, 32);
            for (uint16_t w=0; w < 256; ++w) {
                int8_t excess = 0;
                int8_t rev_excess = 0;
                int32_t min_excess_of_open = 17;
                int32_t min_excess_of_open_pos = 0;
                uint32_t ones = 0;
                min[w] = 8;
                packed_mins[0] = 0x99999999U;
                packed_maxs[0] = 0x99999999U;
                packed_mins.width(4);
                packed_maxs.width(4);
                for (uint16_t p=0; p<8; ++p) {
                    ones += (w&(1<<p))!=0;
                    excess += 1-2*((w&(1<<p))==0);
                    if (excess <= min[w]) {
                        min[w] = excess;
                        min_pos_max[w] = p;
                    }
                    if (excess < 0 and packed_mins[-excess-1] == 9) {
                        packed_mins[-excess-1] = p;
                    }
                    if (w&(1<<p) and excess+8 <= min_excess_of_open) {
                        min_excess_of_open     = excess+8;
                        min_excess_of_open_pos = p;
                    }
                    rev_excess += 1-2*((w&(1<<(7-p)))>0);
                    if (rev_excess < 0 and packed_maxs[-rev_excess-1] == 9) {
                        packed_maxs[-rev_excess-1] = 7-p;
                    }
                }
                word_sum[w] = excess;
                packed_mins.width(32);
                min_match_pos_packed[w] = packed_mins[0];
                packed_maxs.width(32);
                max_match_pos_packed[w] = packed_maxs[0];
                min_open_excess_info[w] = (min_excess_of_open) |
                                          (min_excess_of_open_pos << 8) |
                                          (ones << 12);
            }
        }
    } data;
};

//! Calculate pioneers as defined in the paper of Geary et al. (CPM 2004)
/*! \param bp             The balanced parentheses sequence.
 *  \param block_size     Block size.
 *  \return Bitvector which marks the pioneers in bp.
 *  \par Time complexity
 *       \f$ \Order{n \log n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{2n + min(block\_size, \frac{n}{block\_size} )\cdot \log n } \f$
 */
bit_vector
calculate_pioneers_bitmap(const bit_vector& bp, uint64_t block_size);

//! Space-efficient version of calculate_pioneers_bitmap
/*! \param bp           The balanced parentheses sequence.
 *  \param block_size   Block size.
 *  \return Bitvector which marks the pioneers in bp.
 *  \par Time complexity
 *       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{2n + n} \f$ bits: \f$n\f$ bits for input, \f$n\f$ bits for
 *       output, and \f$n\f$ bits for a succinct stack.
 *  \pre The parentheses sequence represented by bp has to be balanced.
 */
bit_vector
calculate_pioneers_bitmap_succinct(const bit_vector& bp, uint64_t block_size);

//! find_open/find_close for closing/opening parentheses.
/*! \param bp      The balanced parentheses sequence.
 *  \param matches Reference to the result.
 *  \pre bp represents a balanced parentheses sequence.
 *  \par Time complexity
 *       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{n + 2n\log n } \f$
 */
template<class int_vector>
void calculate_matches(const bit_vector& bp, int_vector& matches)
{
    matches = int_vector(bp.size(), 0, bits::hi(bp.size())+1);
    std::stack<uint64_t> opening_parenthesis;
    for (uint64_t i=0; i < bp.size(); ++i) {
        if (bp[i]) {// opening parenthesis
            opening_parenthesis.push(i);
        } else { // closing parenthesis
            assert(!opening_parenthesis.empty());
            uint64_t position = opening_parenthesis.top();
            opening_parenthesis.pop();
            matches[i] = position;
            assert(matches[i]==position);
            matches[position] = i;
            assert(matches[position]==i);
        }
    }
    // assert that the sequence is balanced
    assert(opening_parenthesis.empty());
}

//! Calculates enclose answers for a balanced parentheses sequence.
/*! \param bp A bit_vector representing a balanced parentheses sequence.
 *  \param enclose Reference to the result.
 *  \pre bp represents a balanced parentheses sequence.
 *  \par Time complexity
 *       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
 *  \par Space complexity
 *       \f$ \Order{n + 2n\log n } \f$
 */
template<class int_vector>
void calculate_enclose(const bit_vector& bp, int_vector& enclose)
{
    enclose = int_vector(bp.size(), 0, bits::hi(bp.size())+1);
    std::stack<uint64_t> opening_parenthesis;
    for (uint64_t i=0; i < bp.size(); ++i) {
        if (bp[i]) {// opening parenthesis
            if (!opening_parenthesis.empty()) {
                uint64_t position = opening_parenthesis.top();
                enclose[i] = position;
                assert(enclose[i]==position);
            } else
                enclose[i] = bp.size();
            opening_parenthesis.push(i);
        } else { // closing parenthesis
            uint64_t position = opening_parenthesis.top();
            enclose[i] = position; // find open answer if i is a closing parenthesis
            opening_parenthesis.pop();
        }
    }
    // assert that the sequence is balanced
    assert(opening_parenthesis.empty());
}

uint64_t
near_find_close(const bit_vector& bp, const uint64_t i,
                const uint64_t block_size);

uint64_t
near_find_closing(const bit_vector& bp, uint64_t i,
                  uint64_t closings,
                  const uint64_t block_size);

uint64_t
near_fwd_excess(const bit_vector& bp, uint64_t i,
                bit_vector::difference_type rel, const uint64_t block_size);

//! Calculate the position with minimal excess value in the interval [l..r].
/*! \param bp The bit_vector which represents the parentheses sequence
 *  \param l  The left border of the interval.
 *	\param r  The right border of the interval.
 *  \param min_rel_ex Reference to the relative minimal excess value with regards to excess(bp[l])
 */
uint64_t
near_rmq(const bit_vector& bp, uint64_t l, uint64_t r,
         bit_vector::difference_type& min_rel_ex);

//! Near backward excess
/* This method searches the maximal parenthesis j, with \f$ j\leq i \f$,
 * such that \f$ excess(j) = excess(i+1)+rel \f$ and i < bp.size()-1
 */
uint64_t
near_bwd_excess(const bit_vector& bp, uint64_t i,
                bit_vector::difference_type rel, const uint64_t block_size);

uint64_t
near_find_open(const bit_vector& bp, uint64_t i, const uint64_t block_size);

uint64_t
near_find_opening(const bit_vector& bp, uint64_t i, const uint64_t openings,
                  const uint64_t block_size);

//! Find the opening parenthesis of the enclosing pair if this parenthesis is near.
/*!
 * \param bp bit_vector containing the representation of the balanced parentheses sequence.
 * \param i Position of the opening parenthesis for which we search the position of the opening parenthesis of the enclosing parentheses pair.
 * \param block_size Number of entries to search for the corresponding opening parenthesis of the enclosing parentheses pair.
 * \return If no near enclose exists return i, otherwise the position of the opening parenthesis of the enclosing pair.
 * \pre We assert that \f$ bp[i]=1 \f$
 */
// TODO: implement a fast version using lookup-tables of size 8
uint64_t
near_enclose(const bit_vector& bp, uint64_t i, const uint64_t block_size);

uint64_t
near_rmq_open(const bit_vector& bp, const uint64_t begin, const uint64_t end);

}// end namespace sdsl

#endif
